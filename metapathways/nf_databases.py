"""Database preparation DAG for the public build_db command (no Snakemake)."""
from pathlib import Path
import shlex
import sys
from metapathways.nextflow import task


def metacyc_task(root, source, aligner, memory='16 GB', dependencies=()):
    """Shared by explicit build_db selection and the build_pt dependent export."""
    from metapathways.metacyc_db import TABLES
    root, source = Path(root).expanduser().resolve(), Path(source).expanduser().resolve()
    ext = '.prj' if aligner == 'fast' else '.pdb'
    outputs = [root/'functional/metacyc', root/('functional/formatted/metacyc'+ext),
               root/'functional/formatted/metacyc-names.txt',
               root/'functional_categories/MetaCyc_provenance.json', root/'functional_categories/MetaCyc_reldate.txt']
    outputs.extend(root/'functional_categories'/name for name in TABLES)
    command = shlex.join([sys.executable, '-m', 'metapathways.metacyc_db', '--source', str(source),
                          '--output', str(root), '--aligner', aligner])
    return task('prepare_metacyc', 'Prepare licensed MetaCyc reference', [command], [str(source)],
                list(map(str, outputs)), dependencies, memory=memory, adopt_existing=False,
                cache_outputs=[str(root/'functional/formatted/metacyc.*')])


def plan(root, databases, aligner, test=False, memory='16 GB', metacyc_source=None):
    root = Path(root).resolve()
    scripts = Path(__file__).parent / 'build_DBs'
    q = lambda p: shlex.quote(str(p))
    python = q(sys.executable)
    tasks = []
    dirs = ['functional/formatted', 'functional_categories', 'ncbi_tree', 'taxonomic/formatted', '.metapathways']
    ready = root / '.metapathways/database-directories.ready'
    tasks.append(task('directories', 'Prepare database directories',
        ['mkdir -p ' + ' '.join(q(root/d) for d in dirs) + f'; test -e {q(ready)} || touch {q(ready)}'], outputs=[str(ready)], memory=memory))
    # Recreate missing directories without changing the download dependency timestamp.
    tasks[-1]['status'] = 'redo'

    def add(identifier, commands, inputs, outputs, deps):
        tasks.append(task(identifier, identifier, commands, map(str, inputs), map(str, outputs), deps, memory=memory))
        if identifier.startswith('index_'):
            # Track every index shard so removal of a non-sentinel file invalidates reuse.
            prefix = str(tasks[-1]['outputs'][0]).rsplit('.', 1)[0]
            tasks[-1]['cache_outputs'] = [prefix + '.*']
            tasks[-1]['adopt_existing'] = False
        return identifier

    def fetch(identifier, url, dest, compressed=False):
        dest = Path(dest)
        download = str(dest) + '.download'
        command = f'wget -O {q(download)} {q(url)}\n'
        if compressed:
            command += f'gzip -dc {q(download)} > {q(str(dest)+".tmp")}\nmv {q(str(dest)+".tmp")} {q(dest)}\nrm {q(download)}'
        else:
            command += f'mv {q(download)} {q(dest)}'
        return add(identifier, [command], [ready], [dest], ['directories'])

    if 'metacyc' in databases:
        from metapathways.metacyc_db import source_path
        from metapathways.pt_container import registered_image
        source = metacyc_source or registered_image()
        if not source:
            raise ValueError('MetaCyc requires a licensed source: run build_pt first, or provide --metacyc_source with a complete data directory or SIF. See https://metapathways.readthedocs.io/en/latest/pgdb-workflow.html#metacyc-from-pathway-tools.')
        source = source_path(source)
        tasks.append(metacyc_task(root, source, aligner, memory, ['directories']))
        databases = [db for db in databases if db != 'metacyc']
        if not databases:
            # Adding MetaCyc to an existing MPDB must not refresh unrelated references.
            return tasks

    silvas = ['SILVA_LSU_test', 'SILVA_SSU_test'] if test else [
        'SILVA_138.1_LSURef_NR99_tax_silva_trunc', 'SILVA_138.1_SSURef_NR99_tax_silva_trunc']
    for db in silvas:
        fasta, prefix = root/'taxonomic'/db, root/'taxonomic/formatted'/db
        deps = ['directories']
        if not test:
            deps = [fetch('fetch_'+db, f'https://www.arb-silva.de/fileadmin/silva_databases/release_138_1/Exports/{db}.fasta.gz', fasta, True)]
        add('index_'+db, [f'makeblastdb -in {q(fasta)} -dbtype nucl -parse_seqids -out {q(prefix)}'],
            [fasta], [str(prefix)+'.ndb'], deps)
        add('names_'+db, [f"grep '^>' {q(fasta)} > {q(str(prefix)+'-names.txt')}"], [fasta], [str(prefix)+'-names.txt'], deps)

    if not test:
        note = root/'functional_categories/SILVA_reldate.txt'
        add('release_silva', [f'date -u > {q(note)}\nprintf "Release: 138.1\\n" >> {q(note)}'],
            [], [note], ['fetch_'+db for db in silvas])

    ec = root/'functional_categories/enzyme.dat'
    ec_dep = fetch('fetch_enzyme', 'https://ftp.expasy.org/databases/enzyme/enzyme.dat', ec)
    ec_note = root/'functional_categories/Expasy_reldate.txt'
    add('enzyme_release', [f'wget -O {q(str(ec_note)+".tmp")} https://ftp.expasy.org/databases/enzyme/enzclass.txt\n'
        f'date -u > {q(ec_note)}\ngrep "^Release:" {q(str(ec_note)+".tmp")} >> {q(ec_note)}\nrm {q(str(ec_note)+".tmp")}'],
        [ec], [ec_note], [ec_dep])
    tax_archive = root/'ncbi_tree/new_taxdump.tar.gz'
    tax_dep = fetch('fetch_taxonomy', 'https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz', tax_archive)
    add('taxonomy_tree', [f'tar xzf {q(tax_archive)} -C {q(root/"ncbi_tree")}\n{python} {q(scripts/"taxmap_build.py")} {q(root/"ncbi_tree")}'],
        [tax_archive], [root/'ncbi_tree/ncbi_taxonomy_tree.txt'], [tax_dep])
    unirefs = [db for db in databases if db in ('uniref50', 'uniref90')]
    if unirefs:
        idmap = root/'functional_categories/idmapping.dat.gz'
        idmap_dep = fetch('fetch_idmapping', 'https://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz', idmap)
    for db in databases:
        fasta, prefix = root/'functional'/db, root/'functional/formatted'/db
        deps = ['directories']
        if not test:
            if db == 'swissprot':
                url = 'https://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz'
                deps = [fetch('fetch_'+db, url, fasta, True)]
                fetch('release_'+db, 'https://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/complete/reldate.txt', root/'functional_categories/SwissProt_reldate.txt')
            elif db in unirefs:
                url = f'https://ftp.ebi.ac.uk/pub/databases/uniprot/uniref/{db}/{db}.fasta.gz'
                deps = [fetch('fetch_'+db, url, fasta, True)]
                fetch('release_'+db, f'https://ftp.ebi.ac.uk/pub/databases/uniprot/uniref/{db}/{db}.release_note', root/f'functional_categories/{db}_reldate.txt')
            elif db == 'cazy':
                deps = [fetch('fetch_cazy', 'https://bcb.unl.edu/dbCAN2/download/Databases/V12/CAZyDB.07262023.fa', fasta)]
                note = root/'functional_categories/CAZyDB_reldate.txt'
                add('release_cazy', [f'date -u > {q(note)}\nprintf "Release: V12_07262023\\n" >> {q(note)}'],
                    [fasta], [note], deps)
            elif db == 'eggnog' and not fasta.is_file():
                raise ValueError(f'Provide the eggNOG protein FASTA at {fasta}; the previous database builder did not define an eggNOG download source.')
        ext = '.pdb' if aligner == 'blast' else '.prj'
        cmd = (f'makeblastdb -in {q(fasta)} -dbtype prot -parse_seqids -out {q(prefix)}' if aligner == 'blast'
               else f'fastdb -p {q(prefix)} {q(fasta)}')
        add('index_'+db, [cmd], [fasta], [str(prefix)+ext], deps)
        add('names_'+db, [f"grep '^>' {q(fasta)} > {q(str(prefix)+'-names.txt')}"], [fasta], [str(prefix)+'-names.txt'], deps)
        if db in ('swissprot', 'swissprot_test'):
            add('ec_'+db, [f'{python} {q(scripts/(db+"_mapper.py"))} {q(ec)} {q(root/"functional_categories")}'],
                [fasta, ec], [root/f'functional_categories/EC_map.{db}.tsv'], deps+[ec_dep])
        elif db in unirefs:
            filtered = root/f'functional_categories/{db}_idmapping.dat'
            # A shared download feeds independent maps; no concurrent overwrites.
            add('ec_'+db, [f"gzip -dc {q(idmap)} | grep {q('UniRef'+db[6:]+'_')} > {q(filtered)}\n"
                f'{python} {q(scripts/(db+"_mapper.py"))} {q(ec)} {q(root/"functional_categories")}'],
                [fasta, ec, idmap], [root/f'functional_categories/EC_map.{db}.tsv'], deps+[ec_dep, idmap_dep])
    return tasks

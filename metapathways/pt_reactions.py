"""Filter known unsafe explicit reaction assignments in private PGDB staging."""
import hashlib
import json
from pathlib import Path

BLACKLIST = Path(__file__).parent / 'resources/ptools_reaction_blacklist.json'


def compatibility_path(sample_output):
    if sample_output is None:
        return None
    log = Path(sample_output)/'metapathways_run_log.txt'
    if log.is_file():
        for line in log.read_text().splitlines():
            if line.startswith('Minimum Required Arguments:refdb_dir\t'):
                return Path(line.split('\t', 1)[1])/'functional_categories/ptools_reaction_compatibility.json'
    return None


def image_digest_cached(image, sample_output):
    """Hash once per sample/image identity, even when many MAG workers start together."""
    from metapathways.pt_screen import digest
    if sample_output is None:
        return digest(image)
    import fcntl
    from metapathways.nf_worker import atomic_json
    image = Path(image).resolve()
    stat = image.stat()
    identity = [str(image), stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns]
    cache = Path(sample_output)/'.metapathways/ptools-image-digest.json'
    cache.parent.mkdir(parents=True, exist_ok=True)
    with cache.with_suffix('.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if cache.is_file():
            record = json.loads(cache.read_text())
            if record['identity'] == identity:
                return record['sha256']
        checksum = digest(image)
        atomic_json(cache, dict(identity=identity, sha256=checksum))
        return checksum


def filter_reactions(inputs, image=None, sample_output=None):
    """Keep features and names intact; remove only listed METACYC assignments."""
    inputs = Path(inputs)
    payload = BLACKLIST.read_bytes()
    blacklist = json.loads(payload)
    selected = str(BLACKLIST)
    if image:
        checksum = image_digest_cached(image, sample_output)
        blacklist = {r: entry for r, entry in blacklist.items()
                     if entry.get('image_sha256') == checksum}
    compatibility = compatibility_path(sample_output)
    if compatibility and compatibility.is_file():
        from metapathways.pt_screen import digest
        data = json.loads(compatibility.read_text())
        mapping = compatibility.parent/'MetaCyc-monomer-rxn-pairs.tsv'
        if image and data['image_sha256'] == image_digest_cached(image, sample_output) and mapping.is_file() and data['mapping_sha256'] == digest(mapping):
            payload = compatibility.read_bytes()
            blacklist = data['reactions']
            selected = str(compatibility)
        else:
            print('Pathway Tools compatibility list does not match the image/mapping; using bundled known-trigger fallback', flush=True)
    removals = []
    for file in sorted(inputs.glob('*.pf')):
        feature = None
        kept = []
        changed = False
        for line in file.read_text().splitlines(keepends=True):
            fields = line.rstrip('\r\n').split('\t', 1)
            if fields[0] == 'ID' and len(fields) == 2:
                feature = fields[1]
            if fields[0] == 'METACYC' and len(fields) == 2 and fields[1].strip() in blacklist:
                reaction = fields[1].strip()
                removals.append(dict(file=file.name, feature_id=feature, reaction=reaction,
                                     reason=blacklist[reaction]['reason']))
                changed = True
                continue
            kept.append(line)
            if line.strip() == '//':
                feature = None
        if changed:
            file.write_text(''.join(kept))
    audit = dict(blacklist_source=selected, blacklist_sha256=hashlib.sha256(payload).hexdigest(), removed=removals)
    (inputs/'ptools-reaction-filter.json').write_text(json.dumps(audit, indent=2)+'\n')
    if removals:
        affected = {(r['feature_id'], r['reaction']) for r in removals}
        print(f'Pathway Tools reaction blacklist: removed {len(affected)} explicit feature/reaction '
              'assignments from staged inputs; annotations and feature records retained', flush=True)
    return audit

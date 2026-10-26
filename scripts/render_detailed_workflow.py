"""Render brief and detailed standard flowcharts with the MP palette.

Primary-figure style explicitly requested by the user on 21 October 2026.
The graph describes dependencies, not the order of concurrently scheduled tasks.
"""
from pathlib import Path
import html
import shutil
import re
import cairocffi as cairo
ROOT = Path(__file__).resolve().parents[1]

def render(root, stem=None, output='workflow-detailed', layout_only=False):
    detailed = output == 'workflow-detailed'
    W=1630 if detailed else 1350
    # Compact only the whitespace; shapes, fonts and label padding stay full size.
    def X(x): return 800+(x-900)*.9 if detailed else 660+(x-900)*.74
    def Y(y): return 65+(y-275)*0.85 if detailed else (y if y<=200 else 15+(y-200)*0.65)
    def logical_x(x): return 900+(x-(800 if detailed else 660))/(.9 if detailed else .74)
    metrics=cairo.Context(cairo.ImageSurface(cairo.FORMAT_ARGB32,1,1))
    metrics.select_font_face('Trebuchet MS'); metrics.set_font_size(22)
    p=[]; edges=[]; nodes={}; kinds={}; positions={}
    layout_lines=[]; layout_labels=[]; node_labels={}; input_nodes=set()
    drawing_node=False
    def text(x,y,lines,size=22,mapped=False):
        if not mapped: x,y=X(x),Y(y)
        if not drawing_node: layout_labels.append((x,y,list(lines),size))
        for i,line in enumerate(lines):
            p.append(f'<text x="{x}" y="{y+i*29}" text-anchor="middle" font-size="{size}">{html.escape(line)}</text>')
    def node(key,x,y,lines,kind='process',w=440,h=None,input_db=False):
        nonlocal drawing_node
        lines=[re.sub(r'^[1-4] · ', '', line) for line in lines]
        node_labels[key]=list(lines)
        if kind=='input' or input_db: input_nodes.add(key)
        positions[key]=(x,y)
        x,y=X(x),Y(y)
        # Fit the actual installed font, including padding for slanted sides.
        width=max(metrics.text_extents(line)[2] for line in lines)
        w=max(130,width+(64 if kind in ('input','output','end') else 40))
        h=max(54,len(lines)*29+20)
        if kind=='decision':
            w=2*width+40; h=76
        if kind=='db': h+=24
        nodes[key]=(x,y,w,h); kinds[key]=kind
        fill,stroke=('#DAE8FC','#6C8EBF') if kind=='input' or input_db else ('#D5E8D4','#82B366') if kind in ('output','db','end') else ('#F5F5F5','#666666')
        attrs=f'fill="{fill}" stroke="{stroke}" stroke-width="2" data-node-kind="{kind}"'
        a,b=x-w/2,y-h/2
        if kind in ('input','output'):
            p.append(f'<path d="M{a+20} {b} H{a+w} L{a+w-20} {b+h} H{a} Z" {attrs}/>')
        elif kind=='decision':
            p.append(f'<path d="M{x} {b} L{a+w} {y} L{x} {b+h} L{a} {y} Z" {attrs}/>')
        elif kind=='db':
            # Separate opaque body and lid: a compound path can cancel its
            # own fill where subpaths overlap in SVG/PDF renderers.
            p.append(f'<path d="M{a} {b+12} H{a+w} V{b+h-12} A{w/2} 12 0 0 1 {a} {b+h-12} Z" {attrs}/>')
            p.append(f'<ellipse cx="{x}" cy="{b+12}" rx="{w/2}" ry="12" {attrs}/>')
        else:p.append(f'<rect x="{a}" y="{b}" width="{w}" height="{h}" rx="{h/2 if kind=="end" else 0}" {attrs}/>')
        drawing_node=True
        text(x,y-(len(lines)-1)*14.5+7+(9 if kind=='db' else 0),lines,mapped=True)
        drawing_node=False
    def port(key,side):
        x,y,w,h=nodes[key]
        inset=10 if kinds[key] in ('input','output') else 0
        return {'t':(x,y-h/2),'b':(x,y+h/2),'l':(x-w/2+inset,y),'r':(x+w/2-inset,y),'rlo':(x+w/2,y+min(20,h/2-10))}[side]
    def arrow_path(points,weight=2,dashed=False):
        # End the shaft at the arrow base, not at the tip: SVG markers leave
        # a heavy shaft visible through the taper of a small arrowhead.
        layout_lines.append((list(points),True,dashed,weight))
        pts=list(points)
        ax,ay=pts[-2]; tx,ty=pts[-1]
        length=((tx-ax)**2+(ty-ay)**2)**.5
        if length<=0: raise ValueError('Arrow needs a nonzero final segment')
        ux,uy=(tx-ax)/length,(ty-ay)/length
        head=min(13 if weight>2 else 11,length*.7)
        half=7 if weight>2 else 5
        bx,by=tx-head*ux,ty-head*uy
        pts[-1]=(bx+ux,by+uy)
        edges.append('<path d="'+' '.join(('M' if i==0 else 'L')+f'{x} {y}' for i,(x,y) in enumerate(pts))+f'" fill="none" stroke="#111111" stroke-width="{weight}" stroke-linecap="butt"'+(' stroke-dasharray="7 5"' if dashed else '')+'/>')
        edges.append(f'<polygon points="{tx},{ty} {bx-half*uy},{by+half*ux} {bx+half*uy},{by-half*ux}" fill="#111111" data-arrowhead="true"/>')
    def edge(a,b,sa='b',sb='t',via=(),dashed=False,pixels=False):
        main=((a=='assembly' or a.startswith('s')) and b.startswith('s')) or b=='decision' or (a,b) in {('decision','community'),('extract','pgdb'),('report','end')}
        weight=5.5 if main else 2
        pts=[port(a,sa),*(via if pixels else [(X(x),Y(y)) for x,y in via]),port(b,sb)]
        arrow_path(pts,weight,dashed)
    def rejoin(points,label_x,label_y):
        layout_lines.append((list(points),False,False,2))
        edges.append('<path d="'+' '.join(('M' if i==0 else 'L')+f'{x} {y}' for i,(x,y) in enumerate(points))+'" fill="none" stroke="#111111" stroke-width="2"/>')
        x,y=points[-1]
        p.append(f'<circle cx="{x}" cy="{y}" r="3" fill="#111111"/>')
        text(label_x,label_y,['No'],19,mapped=True)
    # Color encodes role independently of the standard flowchart symbol.
    lx=1060 if detailed else 900
    p.append(f'<rect x="{lx}" y="{20}" width="194" height="{248}" rx="4" fill="#FFFFFF" stroke="#CCCCCC" stroke-width="1.5"/>')
    for i,(label,fill,stroke) in enumerate([('Input','#DAE8FC','#6C8EBF'),('Processing','#F5F5F5','#666666'),('Output','#D5E8D4','#82B366')]):
        cy=40+i*28
        p.append(f'<text x="{lx+114}" y="{cy+7}" text-anchor="end" font-size="20">{label}</text>')
        p.append(f'<rect x="{lx+138}" y="{cy-10}" width="36" height="20" fill="{fill}" stroke="{stroke}" stroke-width="1.5"/>')
    p.append(f'<path d="M{lx+10} 115 H{lx+184}" stroke="#CCCCCC"/>')
    for i,(kind,label) in enumerate([('io','Data I/O'),('process','Process'),('decision','Decision'),('db','Database'),('end','Terminal')]):
        cy=135+i*28; a=lx+128; b=cy-10; w=56
        attrs='fill="#FFFFFF" stroke="#666666" stroke-width="1.5"'
        p.append(f'<text x="{lx+114}" y="{cy+7}" text-anchor="end" font-size="20">{label}</text>')
        if kind=='io':
            p.append(f'<path d="M{a+8} {b} H{a+w} L{a+w-8} {b+20} H{a} Z" {attrs}/>')
        elif kind=='decision':
            p.append(f'<path d="M{a+28} {b-3} L{a+w} {cy} L{a+28} {b+23} L{a} {cy} Z" {attrs}/>')
        elif kind=='db':
            p.append(f'<path d="M{a} {b+4} H{a+w} V{b+17} A28 5 0 0 1 {a} {b+17} Z" {attrs}/>')
            p.append(f'<ellipse cx="{a+28}" cy="{b+4}" rx="28" ry="5" {attrs}/>')
        else: p.append(f'<rect x="{a}" y="{b}" width="{w}" height="20" rx="{10 if kind=="end" else 0}" {attrs}/>')
    node('assembly',900,275,(['Metagenome assemblies','FASTA + sample IDs / manifest'] if detailed else ['Metagenome assemblies']),'input')
    y=430
    stages=[('1 · Preprocessing',[
        ['Validate and filter contigs','Normalize IDs; preserve ID mapping'],
        ['Contig quality control','Sequence statistics; prepared contigs']]),
        ('2 · Feature extraction',[
        ['Predict CDS, rRNA and tRNA','pProdigal · barrnap · ptRNAscan-SE'],
        ['Prepare feature sequences','CDS / proteins and coordinates','Filter protein sequences']]),
        ('3 · Annotation',[
        ['Search protein / rRNA references','FAST / BLAST+; apply hit thresholds'],
        ['Integrate feature annotations','Functional IDs; per-DB taxonomy / LCA','MetaPathways · BEDTools'],
        ['Export annotated features','Tables · GFF / GenBank','Sequence-backed Pathway Tools inputs']])]
    prev='assembly';stage_out={}
    for i,(title,steps) in enumerate(stages,1):
        if not detailed:
            steps=[{1:['Contig QC · MetaPathways'],2:['pProdigal · barrnap · ptRNAscan-SE'],3:['FAST / BLAST+ · taxonomy / LCA']}[i]]
        steps[0] = [title] + steps[0]
        for j,lines in enumerate(steps):
            k=f's{i}_{j}';node(k,900,y,lines);edge(prev,k);prev=k;y+=160
        stage_out[i]=prev
        y+=210 if i==1 else 35
    # Read abundance branches from prepared contigs, then receives feature coordinates.
    prep_y=positions[stage_out[1]][1]
    yscale=.85 if detailed else .65
    mapping_y=prep_y+(nodes[stage_out[1]][3]/2+90)/yscale
    node('reads',260,mapping_y-((78+107)/2+36 if detailed else 90)/yscale,(['Reads','Paired / interleaved FASTQs'] if detailed else ['Reads']),'input',w=420)
    node('mapping',260,mapping_y,['Map reads to prepared contigs','CoverM','Coverage and mapping BAM'] if detailed else ['Map reads · CoverM'],w=420)
    yscale=.85 if detailed else .65
    read_dx=X(900)-175
    read_decision_y=prep_y+(nodes[stage_out[1]][3]/2+90)/yscale
    node('read_decision',logical_x(read_dx),read_decision_y,['Reads?'],'decision')
    edge('reads','mapping')
    split_y=port(stage_out[1],'b')[1]+20
    edge(stage_out[1],'read_decision','b','t',via=[(X(900),split_y),(read_dx,split_y)],pixels=True)
    read_turn=(port('read_decision','l')[0]+port('mapping','r')[0])/2
    edge('read_decision','mapping','l','r')
    text((port('read_decision','l')[0]+port('mapping','r')[0])/2,Y(read_decision_y)-10,['Yes'],19,mapped=True)
    read_join_y=min(port('read_decision','b')[1]+24,port('s2_0','t')[1]-24)
    rejoin([port('read_decision','b'),(read_dx,read_join_y),(X(900),read_join_y)],
           (read_dx+X(900))/2,read_join_y-10)
    anno_y=positions[stage_out[3]][1]
    node('references',1540,positions['s3_0'][1],(['MPDB references','Prepared with build_db'] if detailed else ['MPDB references']),'db',w=420,input_db=True)
    edge('references','s3_0','l','r')
    node('count',260,anno_y,['Count and normalize features','SAMtools · featureCounts · MP','Counts · RPKM · TPM'] if detailed else ['Count / normalize','SAMtools · featureCounts'],w=420)
    edge('mapping','count');edge(stage_out[3],'count','l','r')
    node('abundance',260,anno_y+180,(['Coverage / feature abundance'] if detailed else ['Abundance']),'output',w=420);edge('count','abundance')
    # Optional inference has an explicit bypass; no implied loss of annotations.
    decision_y=anno_y+(nodes[prev][3]/2+170+(107 if detailed else 78)/2+120)/yscale
    node('decision',900,decision_y,['Infer pathways?'],'decision',w=340,h=130);edge(prev,'decision')
    community_x=logical_x(X(900)-(270 if detailed else 230))
    population_x=logical_x(X(900)+(270 if detailed else 230))
    yes_y=port('decision','b')[1]+20
    # Equal whitespace above and below the licensed input.
    license_top= yes_y+32
    license_center=license_top+39
    license_y=decision_y+(license_center-Y(decision_y))/yscale
    node('license',900,license_y,['Licensed','Pathway Tools'],'input')
    build_y=license_y+(39+32+(107 if detailed else 78)/2)/yscale
    # The stage label belongs inside its process box.
    node('community',community_x,build_y,['4 · Community pathway inference','Build PGDB · Pathway Tools / PathoLogic','Compatibility filter; taxonomic pruning'] if detailed else ['4 · Community pathways','Pathway Tools / PathoLogic'])
    yes_y=port('decision','b')[1]+20
    yes_side_x=port('community','l')[0]-28
    edge('decision','community','b','l',via=[(X(900),yes_y),(yes_side_x,yes_y),(yes_side_x,Y(build_y))],pixels=True)
    text((X(900)+yes_side_x)/2,yes_y-10,['Yes'],19,mapped=True)
    edge('license','community','l','t',via=[(X(community_x),Y(license_y))],dashed=True,pixels=True)
    # Genome-map branch reuses community annotations; it does not assemble or bin.
    partition_y=anno_y+(nodes[prev][3]/2+170)/yscale
    node('map',1540,partition_y-(39+36+(107 if detailed else 78)/2)/yscale,['Contig-to-genome','map'],'input',w=420)
    node('partition',1540,partition_y,['Partition annotated features','MAGSplitter','Existing features + genome sequences'] if detailed else ['Partition annotations','MAGSplitter'],w=420)
    map_dx=X(900)+175
    map_decision_y=partition_y
    node('map_decision',logical_x(map_dx),map_decision_y,['MAGs?'],'decision')
    split_y=port(prev,'b')[1]+20
    edge(prev,'map_decision','b','t',via=[(X(900),split_y),(map_dx,split_y)],pixels=True)
    edge('map','partition')
    map_turn=(port('map_decision','r')[0]+min(port('partition','l')[0],port('map','l')[0]-10))/2
    edge('map_decision','partition','r','l')
    text((port('map_decision','r')[0]+port('partition','l')[0])/2,Y(map_decision_y)-10,['Yes'],19,mapped=True)
    map_join_y=max(port(prev,'b')[1],port('map_decision','b')[1])+28
    rejoin([port('map_decision','b'),(map_dx,map_join_y),(X(900),map_join_y)],
           (map_dx+X(900))/2,map_join_y-10)
    # License input sits below partition for both layouts.
    extract_y=build_y+((107 if detailed else 78)/2+64+(78 if detailed else 54)/2)/yscale
    population_y=build_y
    node('population',population_x,population_y,['Build population PGDBs','Pathway Tools / PathoLogic'],w=420)
    edge('partition','population','b','r',via=[(X(1540),Y(population_y))],pixels=True)
    edge('license','population','r','t',via=[(X(population_x),Y(license_y))],dashed=True,pixels=True)
    node('extract',900,extract_y,(['Extract pathways and ORF links','Camelot · MetaPathways'] if detailed else ['Extract pathways · Camelot']))
    def merge_line(points,weight=2):
        layout_lines.append((list(points),False,False,weight))
        edges.append('<path d="'+' '.join(('M' if i==0 else 'L')+f'{x} {y}' for i,(x,y) in enumerate(points))+f'" fill="none" stroke="#111111" stroke-width="{weight}"/>')
    def merged_arrow(y,target):
        arrow_path([(X(900),y),port(target,'t')],5.5)
        p.append(f'<circle cx="{X(900)}" cy="{y}" r="3" fill="#111111"/>')
    extract_join_y=(max(port('community','b')[1],port('population','b')[1])+port('extract','t')[1])/2
    for source in ('community','population'):
        sx,sy=port(source,'b')
        merge_line([(sx,sy),(sx,extract_join_y),(X(900),extract_join_y)],5.5 if source=='community' else 2)
    merged_arrow(extract_join_y,'extract')
    out_y=extract_y+((78 if detailed else 54)/2+36+39)/yscale
    node('pgdb',900,out_y,['Community [population] PGDBs'],'db');edge('extract','pgdb')
    report_y=out_y+(39+64+39)/yscale
    node('report',900,report_y,(['Index all available results','Linked results schema · report / explorer'] if detailed else ['Report / explorer','Linked results schema']))
    join_y=(port('pgdb','b')[1]+port('report','t')[1])/2
    merge_line([port('pgdb','b'),(X(900),join_y)],5.5)
    merged_arrow(join_y,'report')
    edge('abundance','report','b','l',via=[(260,report_y)])
    bypass_x=(X(260)+yes_side_x)/2
    merge_line([port('decision','l'),(bypass_x,Y(decision_y)),(bypass_x,join_y),(X(900),join_y)])
    text(605,decision_y-20,['No'],19)
    node('end',900,report_y+(39+36+39)/yscale,['MP Explorer','Search · subset · export'],'end',w=540);edge('report','end')
    H=port('end','b')[1]+20
    if layout_only:
        return dict(nodes=nodes,kinds=kinds,labels=node_labels,inputs=input_nodes,
                    lines=layout_lines,edge_labels=layout_labels,width=W,height=H)
    header=f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" role="img" aria-labelledby="title desc" data-workflow-style="standard-flowchart"><title id="title">MetaPathways {"detailed" if detailed else "brief"} workflow</title><desc id="desc">Standard flowchart. Rectangles are processes, parallelograms are inputs and outputs, diamonds are decisions, and cylinders are databases. Dashed connectors supply the licensed runtime. Optional population and abundance branches reuse annotated features. Arrows indicate dependencies, not serial execution.</desc><rect width="100%" height="100%" fill="#FFFFFF"/><defs><marker id="arrow" viewBox="0 0 10 10" refX="10" refY="5" markerUnits="userSpaceOnUse" markerWidth="12" markerHeight="12" orient="auto"><path d="M0 0 L10 5 L0 10 Z" fill="#111111"/></marker></defs><g font-family="Trebuchet MS, sans-serif" fill="#111111">'
    target=root/f'docs/assets/{output}.svg';target.write_text(header+'\n'+'\n'.join(edges+p)+'\n</g></svg>\n')
    return target,W,H

def build(root=ROOT):
    for output in ('workflow-main','workflow-detailed'):
        path,_,_=render(root,output=output)
        shutil.copyfile(path,root/'docs/assets/diagrams'/('workflow.svg' if output=='workflow-main' else 'workflow-detailed.svg'))
    return root/'docs/assets/workflow-main.svg'

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pdf',action='store_true')
    args=parser.parse_args();build()
    if args.pdf:
        import cairosvg
        for name in ('workflow-main','workflow-detailed'):
            source=ROOT/'docs/assets'/f'{name}.svg'
            cairosvg.svg2pdf(url=str(source),write_to=str(source.with_suffix('.pdf')))

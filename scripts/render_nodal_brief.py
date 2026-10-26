"""Render the approved MP brief's layout using ASPIRE V1 nodal symbols.

Geometry and connectivity come directly from render_detailed_workflow; this
alternative cannot silently reorder stages or move branches independently.
"""
from pathlib import Path
from html import escape
import math
import cairosvg
from render_detailed_workflow import render, ROOT

TARGET = ROOT / 'docs/assets/workflow-brief-nodal.svg'
# Labels are independent of the shared graph; shared coordinate transforms
# compacts every node and bend together without reorganizing branches.
LABELS = {
 'assembly':(-36,0,'end',['Metagenome assemblies']),
 's1_0':(35,0,'start',['Preprocessing','Contig QC']),
 's2_0':(35,0,'start',['Feature extraction','(Prodigal, barrnap, RNA-scan)']),
 's3_0':(35,0,'start',['Annotation','(FAST, BLAST+, LCA)']),
 'reads':(-35,0,'end',['Reads']),
 'mapping':(-35,0,'end',['Read mapping','(CoverM)']),
 'read_decision':(34,0,'start',['Reads?']),
 'references':(35,0,'start',['MPDB']),
 'count':(-30,0,'end',['Abundance','(featureCounts)']),
 'abundance':(35,0,'start',['Abundance']),
 'decision':(35,0,'start',['Infer pathways?']),
 'license':(0,-46,'middle',['Licensed','Pathway Tools']),
 'community':(28,0,'start',['Community','ePGDB(s)']),
 'map':(35,0,'start',['Contig-to-','genome map']),
 'partition':(35,0,'start',['Partition','(MAGSplitter)']),
 'map_decision':(-34,0,'end',['MAGs?']),
 'population':(-28,0,'end',['Population','ePGDB(s)']),
 'extract':(55,0,'start',['Extract pathways','(Camelot)']),
 'pgdb':(35,0,'start',['Community','[population] ePGDBs']),
 'report':(55,0,'start',['Report / explorer','Linked schema']),
 'end':(35,0,'start',['MP Explorer']),
}

def build():
    layout = render(ROOT, output='workflow-main', layout_only=True)
    # Feed references into the spine between feature extraction and annotation,
    # keeping the annotation label clear of the horizontal connector.
    rx, ry, rw, rh = layout['nodes']['references']
    annotation_x, annotation_y, _, _ = layout['nodes']['s3_0']
    join_y = (layout['nodes']['s2_0'][1] + annotation_y) / 2
    layout['nodes']['references'] = (rx, join_y, rw, rh)
    layout['lines'] = [
        ([(rx-rw/2, join_y), (annotation_x, join_y)], arrow, dashed, weight)
        if points[0] == (rx-rw/2, ry) else (points, arrow, dashed, weight)
        for points, arrow, dashed, weight in layout['lines']
    ]
    # Lower the license input and its supply rails together.
    lx, ly, lw, lh = layout['nodes']['license']
    license_drop = 37 / .65
    layout['nodes']['license'] = (lx, ly+license_drop, lw, lh)
    layout['lines'] = [
        ([(x, y+license_drop if abs(y-ly)<1e-8 else y) for x,y in points],
         arrow, dashed, weight) if dashed else (points, arrow, dashed, weight)
        for points, arrow, dashed, weight in layout['lines']
    ]
    edges, shapes = [], []
    def compact_x(x):
        # Pull the outer lanes inward without narrowing the paired ePGDB builds.
        shift = 65 * min(1, max(0, (430-x)/(430-244)))
        shift -= 60 * min(1, max(0, (x-890)/(1133.6-890)))
        return 180+(x-186.4)*.60+shift-60
    # Preserve the supply-arrow clearance after lowering the license input.
    vertical_padding=[
        # Open equal space above and below MPDB without shifting feature text.
        ((layout['nodes']['s2_0'][1]+join_y)/2,18),
        ((join_y+annotation_y)/2,18),
        (layout['nodes']['license'][1]+20,12),
    ]
    def compact_y(y):
        return 40+(y-63.75)*.65+sum(gap for threshold,gap in vertical_padding if y>=threshold)
    def label(x,y,lines,anchor='middle',size=20,software=False):
        for i,value in enumerate(lines):
            yy=y+7+(i-(len(lines)-1)/2)*25
            line_size=17 if software and value.startswith('(') else size
            shapes.append(f'<text x="{x}" y="{yy}" font-size="{line_size}" text-anchor="{anchor}">{escape(value)}</text>')
    old_ports=set()
    for key,(x,y,w,h) in layout['nodes'].items():
        inset=10 if layout['kinds'][key] in ('input','output') else 0
        old_ports.update(((x,y-h/2),(x,y+h/2),(x-w/2+inset,y),
                          (x+w/2-inset,y),(x+w/2,y+min(20,h/2-10))))
    # A junction is an endpoint/bend with at least three incident directions.
    # Mark genuine joins and splits, not ordinary bends or unconnected crossings.
    segments = [(a,b) for points,_,_,_ in layout['lines']
                for a,b in zip(points,points[1:])]
    candidates = {point for points,_,_,_ in layout['lines'] for point in points}
    junctions = set()
    for point in candidates - old_ports:
        px,py = point
        directions = set()
        for (ax,ay),(bx,by) in segments:
            cross = (px-ax)*(by-ay)-(py-ay)*(bx-ax)
            if abs(cross)>1e-7 or not (min(ax,bx)-1e-7<=px<=max(ax,bx)+1e-7
                                      and min(ay,by)-1e-7<=py<=max(ay,by)+1e-7):
                continue
            for x,y in ((ax,ay),(bx,by)):
                length = math.hypot(x-px,y-py)
                if length>1e-7:
                    directions.add((round((x-px)/length,6),round((y-py)/length,6)))
        if len(directions)>=3:
            junctions.add(point)
    # Preserve a visible shaft between junction dots and downstream heads.
    # Add space only in short affected intervals, retaining the compact layout.
    for key,(x,y,w,h) in sorted(layout['nodes'].items(),key=lambda item:item[1][1]):
        radius=18 if layout['kinds'][key] in ('input','output','db','end') else 21
        for points,arrow,_,weight in layout['lines']:
            if not arrow or points[-1]!=(x,y-h/2) or points[-2][0]!=x:
                continue
            upstream=[jy for jx,jy in junctions if jx==x and points[-2][1]<=jy<y-h/2]
            if not upstream:
                continue
            jy=max(upstream)
            required=radius+(13 if weight>2 else 9)+4.5+12
            extra=required-(compact_y(y)-compact_y(jy))
            if extra>0:
                vertical_padding.append(((jy+y)/2,extra))
    port_map={}
    for key,(x,y,w,h) in layout['nodes'].items():
        kind=layout['kinds'][key]
        radius=18 if kind in ('input','output','db','end') else 21
        inset=10 if kind in ('input','output') else 0
        for old,new in [((x,y-h/2),(x,y-radius)),((x,y+h/2),(x,y+radius)),
                        ((x-w/2+inset,y),(x-radius,y)),((x+w/2-inset,y),(x+radius,y)),
                        ((x+w/2,y+min(20,h/2-10)),(x+radius,y))]:
            port_map[old]=(compact_x(x)+new[0]-x,compact_y(y)+new[1]-y)
        x,y=compact_x(x),compact_y(y)
        shapes.append(f'<g id="node-{key}" data-center-x="{x}" data-center-y="{y}">')
        if radius==18:
            fill,stroke=('#DAE8FC','#6C8EBF') if key in layout['inputs'] else ('#D5E8D4','#82B366')
            shapes.append(f'<circle cx="{x}" cy="{y}" r="18" fill="{fill}" stroke="{stroke}" stroke-width="3"/>')
        else:
            shapes.append(f'<path d="M{x} {y-21} L{x+21} {y} L{x} {y+21} L{x-21} {y} Z" fill="white" stroke="#111111" stroke-width="1.2"/>')
        dx,dy,anchor,lines=LABELS[key]
        label(x+dx,y+dy,lines,anchor,software=key!='license')
        shapes.append('</g>')
    for original,arrow,dashed,weight in layout['lines']:
        arrow = arrow and original[-1] not in junctions
        width=weight if weight>2 else 1.4
        pts=list(original)
        # Only ports move to the smaller symbols. Every branch bend is retained.
        pts=[(compact_x(x),compact_y(y)) for x,y in pts]
        pts[0]=port_map.get(original[0],pts[0]); pts[-1]=port_map.get(original[-1],pts[-1])
        if arrow:
            ax,ay=pts[-2]; x,y=pts[-1]
            length=math.hypot(x-ax,y-ay)
            ux,uy=(x-ax)/length,(y-ay)/length
            head=min(13 if weight>2 else 9,length*.7); bx,by=x-head*ux,y-head*uy
            half=7 if weight>2 else 4
            pts[-1]=(bx,by)
            edges.append(f'<path d="M{x} {y} L{bx-half*uy} {by+half*ux} L{bx+half*uy} {by-half*ux} Z" fill="#111111"/>')
        d=' '.join(('M' if i==0 else 'L')+f'{x} {y}' for i,(x,y) in enumerate(pts))
        edges.append(f'<path d="{d}" fill="none" stroke="#111111" stroke-width="{width}" stroke-linecap="butt"'+(' stroke-dasharray="6 5"' if dashed else '')+'/>')
    for x,y in sorted(junctions):
        shapes.append(f'<circle cx="{compact_x(x)}" cy="{compact_y(y)}" r="4.5" fill="#111111" data-junction="true"/>')
    for x,y,lines,size in layout['edge_labels']:
        # The two optional-branch No labels belong below their return rails.
        below = lines == ['No'] and x > 500
        label(compact_x(x),compact_y(y+10)+17 if below else compact_y(y)-7,
              lines,size=size)
    # Same legend location; the nodal key uses the ASPIRE V1 three symbols.
    for y,role in [(35,'Input'),(77,'Process'),(119,'Output')]:
        x=720
        if role=='Process':
            shapes.append(f'<path d="M{x} {y-17} L{x+17} {y} L{x} {y+17} L{x-17} {y} Z" fill="white" stroke="#111111"/>')
        else:
            fill,stroke=('#DAE8FC','#6C8EBF') if role=='Input' else ('#D5E8D4','#82B366')
            shapes.append(f'<circle cx="{x}" cy="{y}" r="15" fill="{fill}" stroke="{stroke}" stroke-width="3"/>')
        label(x-35,y,[role],'end')
    w,h=820,compact_y(layout['height'])
    svg=f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}" role="img" data-workflow-style="aspire-v1-nodal"><title>MetaPathways brief workflow</title><desc>Branch organization and connector bends from the approved brief, compressed horizontally and vertically and rendered with ASPIRE V1-style symbols and external labels.</desc><rect width="100%" height="100%" fill="white"/><g font-family="Trebuchet MS, sans-serif" fill="#111111">'+'\n'.join(edges+shapes)+'</g></svg>'
    TARGET.write_text(svg)
    return TARGET

if __name__=='__main__':
    source=build()
    cairosvg.svg2pdf(url=str(source),write_to=str(source.with_suffix('.pdf')))

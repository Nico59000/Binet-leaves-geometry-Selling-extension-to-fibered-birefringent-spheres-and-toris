#!/usr/bin/env python3
from pathlib import Path
import fitz, hashlib, json, math
from PIL import Image, ImageOps, ImageDraw
import numpy as np
B=Path(__file__).resolve().parent
pre=B/'v53_precanon.pdf'; can=B/'formalisation_feuillets_binet_birefringence_selling_polaire_v53.pdf'
d0=fitz.open(pre); d1=fitz.open(can)
assert len(d0)==len(d1)
DPI=60; mat=fitz.Matrix(DPI/72,DPI/72)
diff=[]; flags=[]; min_ink=1.; max_ink=0.; max_dark=0.; max_border=0.
heads={'v53':[],'certification':[],'conclusion':[],'provenance':[],'bibliography':[]}
keywords={'v53':'v53 :','certification':'Certification HT+NT','conclusion':'Conclusion','provenance':'Provenance','bibliography':'Références documentaires du corpus'}
for i,(p0,p1) in enumerate(zip(d0,d1),start=1):
    a=np.frombuffer(p0.get_pixmap(matrix=mat,alpha=False,colorspace=fitz.csGRAY).samples,dtype=np.uint8)
    b=np.frombuffer(p1.get_pixmap(matrix=mat,alpha=False,colorspace=fitz.csGRAY).samples,dtype=np.uint8)
    if a.shape!=b.shape or not np.array_equal(a,b): diff.append(i)
    arr=b.reshape(p1.get_pixmap(matrix=mat,alpha=False,colorspace=fitz.csGRAY).height,-1)
    ink=float((arr<245).mean()); dark=float((arr<80).mean())
    bw=max(1,min(arr.shape)//100); border=np.concatenate([arr[:bw,:].ravel(),arr[-bw:,:].ravel(),arr[:, :bw].ravel(),arr[:, -bw:].ravel()])
    bd=float((border<80).mean())
    min_ink=min(min_ink,ink);max_ink=max(max_ink,ink);max_dark=max(max_dark,dark);max_border=max(max_border,bd)
    if ink<0.01 or ink>0.35 or dark>0.08 or bd>0.03: flags.append({'page':i,'ink':ink,'dark':dark,'border_dark':bd})
    t=p1.get_text()
    for k,kw in keywords.items():
        if kw in t: heads[k].append(i)
# contact sheets around v53 and terminal pages
ranges=[(260,268),(269,len(d1))]
contacts=[]
for lo,hi in ranges:
    thumbs=[]
    for pn in range(lo,hi+1):
        pix=d1[pn-1].get_pixmap(matrix=fitz.Matrix(0.8,0.8),alpha=False)
        im=Image.frombytes('RGB',[pix.width,pix.height],pix.samples)
        im.thumbnail((420,594))
        canvas=Image.new('RGB',(440,625),'white');canvas.paste(im,((440-im.width)//2,20));
        ImageDraw.Draw(canvas).text((10,5),f'p.{pn}',fill='black');thumbs.append(canvas)
    cols=3;rows=math.ceil(len(thumbs)/cols);sheet=Image.new('RGB',(cols*440,rows*625),'white')
    for j,im in enumerate(thumbs):sheet.paste(im,((j%cols)*440,(j//cols)*625))
    fn=B/f'contact_v53_{lo}_{hi}.jpg';sheet.save(fn,quality=88);contacts.append(fn.name)
report={'version':'v53','status':'PASS' if (not diff and not flags) else 'FAIL','pdf_pages':len(d1),'all_pages_rendered':True,
'precanon_canonical_all_page_pixel_parity':{'pages':len(d1),'identical_pages':len(d1)-len(diff),'different_pages':diff,'status':'PASS' if not diff else 'FAIL'},
'automated_render_scan':{'pages':len(d1),'flags':flags,'min_ink':min_ink,'max_ink':max_ink,'max_dark':max_dark,'max_border_dark':max_border,'status':'PASS' if not flags else 'FAIL'},
'pdf_pages_by_heading':heads,'visual_contact_sheets':contacts,'renderer':'PyMuPDF grayscale at 60 dpi for full-page parity and scan; contact sheets at 57.6 dpi'}
(B/'v53_PDF_QA_RAW.json').write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n')
print(json.dumps(report,indent=2,ensure_ascii=False))

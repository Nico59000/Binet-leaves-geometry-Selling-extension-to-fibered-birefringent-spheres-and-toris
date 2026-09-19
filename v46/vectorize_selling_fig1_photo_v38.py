#!/usr/bin/env python3
from pathlib import Path
from collections import defaultdict, Counter
import hashlib, json, math
import cv2, numpy as np
from skimage.morphology import skeletonize
from scipy.ndimage import convolve
from scipy.spatial import cKDTree

ROOT=Path('/mnt/data/formalisation_reprise')
SRC=ROOT/'Tafel-III-fig1.jpg'
OLD=ROOT/'selling_fig1_full_600.png'
OLD_GEO=ROOT/'selling_fig1_layers_v37.geojson'
OUT=ROOT
CROP_BOX=(1460,120,3550,3000) # x0,y0,x1,y1 in full-photo pixels

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def sample_polyline(pts,step=2.0):
    pts=np.asarray(pts,float)
    out=[]
    for a,b in zip(pts[:-1],pts[1:]):
        L=float(np.linalg.norm(b-a)); n=max(2,int(math.ceil(L/step))+1)
        for t in np.linspace(0,1,n,endpoint=False): out.append(a+t*(b-a))
    if len(pts): out.append(pts[-1])
    return np.asarray(out,float)

def transform_points(H,pts):
    P=np.asarray(pts,float)
    ph=np.c_[P,np.ones(len(P))]
    q=(H@ph.T).T
    return q[:,:2]/q[:,2,None]

full=cv2.imread(str(SRC),cv2.IMREAD_COLOR)
assert full is not None
x0,y0,x1,y1=CROP_BOX
crop=full[y0:y1,x0:x1].copy()
crop_gray=cv2.cvtColor(crop,cv2.COLOR_BGR2GRAY)
cv2.imwrite(str(OUT/'selling_fig1_photo_crop_pre_v38.png'),crop)

old600=cv2.imread(str(OLD),cv2.IMREAD_GRAYSCALE)
assert old600 is not None and old600.shape==(5400,4450)
old=cv2.resize(old600,(2225,2700),interpolation=cv2.INTER_AREA)

# Source registration: photo crop -> v37 300-dpi work image.
def enhance(im):
    return cv2.createCLAHE(clipLimit=2.0,tileGridSize=(16,16)).apply(im)
sift=cv2.SIFT_create(nfeatures=8000,contrastThreshold=0.015,edgeThreshold=20)
kp_n,des_n=sift.detectAndCompute(enhance(crop_gray),None)
kp_o,des_o=sift.detectAndCompute(enhance(old),None)
bf=cv2.BFMatcher(cv2.NORM_L2)
knn=bf.knnMatch(des_n,des_o,k=2)
good=[m for m,n in knn if m.distance < 0.70*n.distance]
src=np.float32([kp_n[m.queryIdx].pt for m in good]).reshape(-1,1,2)
dst=np.float32([kp_o[m.trainIdx].pt for m in good]).reshape(-1,1,2)
H,inl=cv2.findHomography(src,dst,cv2.RANSAC,4.0)
assert H is not None
inmask=inl.ravel().astype(bool)
src_in=src.reshape(-1,2)[inmask]; dst_in=dst.reshape(-1,2)[inmask]
proj=transform_points(H,src_in)
err=np.linalg.norm(proj-dst_in,axis=1)
Hinv=np.linalg.inv(H)

# Background-normalized darkness from the photo; this retains weak grey strokes.
bg=cv2.GaussianBlur(crop_gray,(0,0),sigmaX=35,sigmaY=35)
norm=cv2.divide(crop_gray,bg,scale=255)
b=(norm<235).astype(np.uint8)
num,lab,stats,cent=cv2.connectedComponentsWithStats(b,8)
keep=np.zeros_like(b)
for i in range(1,num):
    x,y,w,h,a=map(int,stats[i]); mx=max(w,h); mn=min(w,h)
    elongated=mx>=55 and a>=25 and (mx/(mn+1)>=1.6 or a>=180)
    if elongated: keep[lab==i]=1
# Source-specific nuisance masks: page gutter, figure title, lithographer signature.
keep[:,:120]=0
keep[:520,:280]=0
keep[2500:,1450:]=0
keep=cv2.morphologyEx(keep,cv2.MORPH_CLOSE,np.ones((3,3),np.uint8),iterations=1)
cv2.imwrite(str(OUT/'selling_fig1_photo_longstroke_mask_pre_v38.png'),np.where(keep,0,255).astype(np.uint8))

# Skeleton graph, intentionally diagnostic rather than semantic.
sk=skeletonize(keep.astype(bool))
ker=np.ones((3,3),np.uint8); ker[1,1]=0
degpix=convolve(sk.astype(np.uint8),ker,mode='constant',cval=0)
nodepix=sk & (degpix!=2)
npix=cv2.dilate(nodepix.astype(np.uint8),np.ones((3,3),np.uint8),iterations=1)
nn,nlab,nstats,ncent=cv2.connectedComponentsWithStats(npix,8)
node_id=np.zeros_like(nlab,np.int32)
for j,i in enumerate(range(1,nn),start=1): node_id[nlab==i]=j
edgepix=sk & (node_id==0)
ne,elab,estats,ecent=cv2.connectedComponentsWithStats(edgepix.astype(np.uint8),8)
nbr=[(-1,-1),(-1,0),(-1,1),(0,-1),(0,1),(1,-1),(1,0),(1,1)]
def order_component(coords):
    S=set(map(tuple,coords.tolist()))
    adj={p:[(p[0]+dy,p[1]+dx) for dy,dx in nbr if (p[0]+dy,p[1]+dx) in S] for p in S}
    ends=[p for p,v in adj.items() if len(v)<=1]
    start=ends[0] if ends else next(iter(S)); path=[start]; prev=None; cur=start; seen={cur}
    for _ in range(len(S)+5):
        cand=[q for q in adj[cur] if q!=prev and q not in seen]
        if not cand: break
        q=cand[0]; path.append(q); seen.add(q); prev,cur=cur,q
    return path
DT=cv2.distanceTransform(b,cv2.DIST_L2,5)
raw=[]
for i in range(1,ne):
    x,y,w,h,a=map(int,estats[i])
    if a<8: continue
    ys,xs=np.where(elab==i); coords=np.column_stack([ys,xs])
    comp=(elab==i).astype(np.uint8)
    ring=cv2.dilate(comp,np.ones((3,3),np.uint8),1).astype(bool) & ~comp.astype(bool)
    ids=sorted(set(int(v) for v in node_id[ring] if v>0))
    if len(ids)!=2: continue
    path=order_component(coords)
    xy=np.array([[p[1],p[0]] for p in path],dtype=np.float32)
    if len(xy)<2: continue
    plen=float(np.sum(np.hypot(np.diff(xy[:,0]),np.diff(xy[:,1]))))
    if plen<35: continue
    eps=max(1.2,min(4.0,plen/100))
    ap=cv2.approxPolyDP(xy.reshape(-1,1,2),eps,False).reshape(-1,2)
    chord=float(np.linalg.norm(ap[-1]-ap[0])) if len(ap)>1 else 0.0
    tort=plen/max(chord,1e-9)
    samp=[]
    for px,py in xy[::max(1,len(xy)//150)]:
        xx=int(round(float(px))); yy=int(round(float(py)))
        if 2<=xx<crop_gray.shape[1]-2 and 2<=yy<crop_gray.shape[0]-2:
            samp.append(float(np.max(DT[yy-1:yy+2,xx-1:xx+2])))
    raw.append({'raw_edge':i-1,'node_ids':[j-1 for j in ids], 'length_px_photo':plen,
                'tortuosity':tort,'stroke_radius_med_px_photo':float(np.median(samp)) if samp else 0.0,
                'polyline':ap.astype(float)})

oldgeo=json.loads(OLD_GEO.read_text(encoding='utf-8'))
old_samples=[]
for f in oldgeo['features']:
    if f['properties']['layer']=='symmetry_candidate': continue
    q=transform_points(Hinv,f['geometry']['coordinates'])
    old_samples.append(sample_polyline(q,1.5))
tree=cKDTree(np.vstack(old_samples))

features=[]; counts=Counter()
for j,e in enumerate(raw):
    sm=sample_polyline(e['polyline'],2.0)
    d,_=tree.query(sm,k=1)
    f6=float(np.mean(d<=6)); f10=float(np.mean(d<=10))
    if f6>0.50:
        cls='PHOTO_MATCHED_TO_V37_GEOMETRY'
    elif f6<0.20 and e['length_px_photo']>=70:
        cls='PHOTO_NEW_LONG_GEOMETRY_CANDIDATE'
    else:
        cls='PHOTO_PARTIAL_OR_FRAGMENT_CANDIDATE'
    if e['stroke_radius_med_px_photo']>=2.8 and e['length_px_photo']>=50:
        style='REINFORCED_STROKE_CANDIDATE'
    elif e['tortuosity']>=1.18 and e['length_px_photo']>=60:
        style='CURVED_OR_SERPENTINE_CANDIDATE'
    else:
        style='ORDINARY_STROKE_CANDIDATE'
    counts[cls]+=1; counts[style]+=1
    props={'candidate_id':f'P{j:03d}','class':cls,'style_candidate':style,
           'length_px_photo':round(e['length_px_photo'],3),'tortuosity':round(e['tortuosity'],6),
           'stroke_radius_med_px_photo':round(e['stroke_radius_med_px_photo'],3),
           'v37_support_fraction_6px':round(f6,6),'v37_support_fraction_10px':round(f10,6),
           'status':'DIAGNOSTIC_ONLY_NOT_HISTORICAL_INCIDENCE'}
    features.append({'type':'Feature','geometry':{'type':'LineString','coordinates':[[round(float(x),2),round(float(y),2)] for x,y in e['polyline']]},'properties':props})

geo_out={'type':'FeatureCollection','crs_note':'pixel coordinates in canonical photo crop; origin top-left',
         'source':SRC.name,'source_sha256':sha(SRC),'crop_box_full_photo':CROP_BOX,
         'status':'PRE_V38_DIAGNOSTIC_ONLY','features':features}
(OUT/'selling_fig1_photo_candidates_pre_v38.geojson').write_text(json.dumps(geo_out,indent=2,ensure_ascii=False),encoding='utf-8')

# Visual cross-check: v37 geometry transferred into photo coordinates plus photo candidates.
ov=crop.copy()
# transferred v37 source strokes
for f in oldgeo['features']:
    q=transform_points(Hinv,f['geometry']['coordinates'])
    pts=np.round(q).astype(np.int32).reshape(-1,1,2)
    if f['properties']['layer']=='symmetry_candidate': col=(200,120,0)
    else: col=(220,40,40)
    cv2.polylines(ov,[pts],False,col,1,cv2.LINE_AA)
for feat in features:
    pts=np.round(np.array(feat['geometry']['coordinates'])).astype(np.int32).reshape(-1,1,2)
    cls=feat['properties']['class']
    if cls=='PHOTO_MATCHED_TO_V37_GEOMETRY': col=(20,150,20)
    elif cls=='PHOTO_NEW_LONG_GEOMETRY_CANDIDATE': col=(20,20,220)
    else: col=(0,140,220)
    cv2.polylines(ov,[pts],False,col,1,cv2.LINE_AA)
cv2.imwrite(str(OUT/'selling_fig1_photo_crosscheck_pre_v38.png'),ov)

# Support of each v37 layer in the new photo.
bg2=cv2.GaussianBlur(crop_gray,(0,0),sigmaX=18,sigmaY=18)
dark=np.clip(bg2.astype(np.int16)-crop_gray.astype(np.int16),0,255).astype(np.uint8)
inkd=cv2.dilate((dark>12).astype(np.uint8),np.ones((5,5),np.uint8),1)
support=defaultdict(list)
for f in oldgeo['features']:
    q=transform_points(Hinv,f['geometry']['coordinates']); sm=sample_polyline(q,2.0); vals=[]
    for x,y in sm:
        xx=int(round(x)); yy=int(round(y))
        if 0<=xx<inkd.shape[1] and 0<=yy<inkd.shape[0]: vals.append(bool(inkd[yy,xx]))
    support[f['properties']['layer']].append(float(np.mean(vals)) if vals else 0.0)
support_summary={k:{'n':len(v),'median':float(np.median(v)),'mean':float(np.mean(v)),'min':float(np.min(v)),'max':float(np.max(v))} for k,v in support.items()}

manifest={
  'status':'PRE_V38_SOURCE_REGISTRATION_AND_PHOTO_VECTOR_DIAGNOSTIC_NOT_A_CYCLE_CLOSE',
  'new_source':{'file':SRC.name,'sha256':sha(SRC),'bytes':SRC.stat().st_size,'dimensions_full':[int(full.shape[1]),int(full.shape[0])],
                'canonical_crop_box_xyxy':list(CROP_BOX),'canonical_crop_file':'selling_fig1_photo_crop_pre_v38.png',
                'canonical_crop_sha256':sha(OUT/'selling_fig1_photo_crop_pre_v38.png'),'canonical_crop_dimensions':[int(crop.shape[1]),int(crop.shape[0])]},
  'v37_reference':{'file':OLD.name,'sha256':sha(OLD),'working_dimensions_300':[2225,2700]},
  'registration':{'method':'SIFT + Lowe 0.70 + RANSAC homography','raw_good_matches':len(good),'inliers':int(inmask.sum()),
                  'inlier_fraction':float(inmask.mean()),'median_reprojection_error_px_v37_300':float(np.median(err)),
                  'mean_reprojection_error_px_v37_300':float(np.mean(err)),'p95_reprojection_error_px_v37_300':float(np.percentile(err,95)),
                  'H_photo_to_v37_300':H.tolist()},
  'photo_vector_diagnostic':{'threshold_normalized_lt':235,'candidate_edge_count':len(features),'class_counts':dict(counts),
                             'guard':'Counts are algorithmic stroke fragments/candidates, not the historical field/segment count.'},
  'v37_geometry_support_in_photo':support_summary,
  'decision':'PHOTO_SUPERSEDES_RASTER_FOR_GEOMETRIC_EXTRACTION; V37_RETAINED_AS_CROSSCHECK; HISTORICAL_INCIDENCE_REMAINS_NT_UNTIL_SEMANTIC_MATCHING'
}
(OUT/'selling_fig1_photo_registration_pre_v38.json').write_text(json.dumps(manifest,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps(manifest,indent=2,ensure_ascii=False))

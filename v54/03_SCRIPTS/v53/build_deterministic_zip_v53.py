#!/usr/bin/env python3
from pathlib import Path
import zipfile,sys
root=Path(sys.argv[1]).resolve(); out=Path(sys.argv[2]).resolve()
files=sorted(p for p in root.rglob('*') if p.is_file() and p.resolve()!=out)
with zipfile.ZipFile(out,'w',compression=zipfile.ZIP_STORED,allowZip64=True) as z:
    for p in files:
        rel=p.relative_to(root).as_posix()
        zi=zipfile.ZipInfo(rel,date_time=(1980,1,1,0,0,0));zi.compress_type=zipfile.ZIP_STORED
        zi.external_attr=(0o100644 & 0xFFFF)<<16;zi.create_system=3
        z.writestr(zi,p.read_bytes(),compress_type=zipfile.ZIP_STORED)
print(len(files))

#!/usr/bin/env python3
from pypdf import PdfWriter
from pypdf.generic import ArrayObject, ByteStringObject
from pathlib import Path
import sys,tempfile,os
src=Path(sys.argv[1]);out=Path(sys.argv[2]);fixed=bytes.fromhex('76353300000000000000000000000000')
def once(inp,op):
    w=PdfWriter(clone_from=str(inp));w.metadata=None;w._ID=ArrayObject([ByteStringObject(fixed),ByteStringObject(fixed)])
    with open(op,'wb') as f:w.write(f)
fd,tmp=tempfile.mkstemp(prefix='v53canon_',suffix='.pdf',dir=str(out.parent));os.close(fd)
try:
    once(src,tmp);once(tmp,out)
finally:
    try:os.remove(tmp)
    except FileNotFoundError:pass

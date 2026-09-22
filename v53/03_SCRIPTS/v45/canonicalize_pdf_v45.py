#!/usr/bin/env python3
from pypdf import PdfWriter
from pypdf.generic import ArrayObject, ByteStringObject
from pathlib import Path
import sys
src=Path(sys.argv[1]); out=Path(sys.argv[2])
w=PdfWriter(clone_from=str(src)); w.metadata=None
fixed=bytes.fromhex('76343500000000000000000000000000')
w._ID=ArrayObject([ByteStringObject(fixed),ByteStringObject(fixed)])
with out.open('wb') as f:w.write(f)

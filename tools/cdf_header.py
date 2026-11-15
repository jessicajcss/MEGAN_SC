"""Minimal parser for the NetCDF classic / 64-bit-offset / CDF-5 *header*.

Reads dims, variables (name/type/shape/dims) and attributes from just the front of a
file — so it works on a truncated slice of a multi-GB wrfout without netCDF4. Spec:
https://docs.unidata.ucar.edu/netcdf-c/current/file_format_specifications.html
"""
from __future__ import annotations

import struct

NC_TYPE = {1: "byte", 2: "char", 3: "short", 4: "int", 5: "float", 6: "double",
           7: "ubyte", 8: "ushort", 9: "uint", 10: "int64", 11: "uint64"}
_TSIZE = {1: 1, 2: 1, 3: 2, 4: 4, 5: 4, 6: 8, 7: 1, 8: 2, 9: 4, 10: 8, 11: 8}


class R:
    def __init__(self, b): self.b, self.p = b, 0
    def u32(self):
        v = struct.unpack_from(">I", self.b, self.p)[0]; self.p += 4; return v
    def u64(self):
        v = struct.unpack_from(">Q", self.b, self.p)[0]; self.p += 8; return v
    def name(self):
        n = self.u32(); s = self.b[self.p:self.p + n].decode("latin-1", "replace")
        self.p += n + (-n % 4); return s
    def values(self, nc_type, nelems):
        sz = _TSIZE[nc_type]; raw = self.b[self.p:self.p + sz * nelems]
        self.p += sz * nelems + (-(sz * nelems) % 4)
        if nc_type == 2:
            return raw.decode("latin-1", "replace").rstrip("\x00 ")
        fmt = {1: "b", 3: "h", 4: "i", 5: "f", 6: "d", 7: "B", 8: "H", 9: "I",
               10: "q", 11: "Q"}[nc_type]
        return list(struct.unpack(">" + fmt * nelems, raw))


def parse(path):
    with open(path, "rb") as f:
        b = f.read()
    magic, ver = b[:3], b[3]
    assert magic == b"CDF", f"not a classic NetCDF: {b[:4]!r}"
    off64 = ver == 2
    r = R(b); r.p = 4
    numrecs = r.u32() if ver != 5 else r.u64()

    def dimsize():   return r.u64() if ver == 5 else r.u32()
    def offset():    return r.u64() if (off64 or ver == 5) else r.u32()

    def att_list():
        tag = r.u32(); n = r.u32() if ver != 5 else r.u64()
        atts = {}
        if tag == 0 and n == 0:
            return atts
        for _ in range(n):
            nm = r.name(); tp = r.u32(); ne = (r.u64() if ver == 5 else r.u32())
            atts[nm] = r.values(tp, ne)
        return atts

    # dimensions
    tag = r.u32(); nd = r.u32() if ver != 5 else r.u64()
    dims = []
    for _ in range(nd if tag else 0):
        dims.append((r.name(), dimsize()))
    gatts = att_list()

    # variables
    tag = r.u32(); nv = r.u32() if ver != 5 else r.u64()
    vars_ = []
    for _ in range(nv if tag else 0):
        nm = r.name(); ndims = r.u32() if ver != 5 else r.u64()
        dimids = [(r.u64() if ver == 5 else r.u32()) for _ in range(ndims)]
        vatts = att_list(); tp = r.u32()
        vsize = r.u64() if ver == 5 else r.u32(); begin = offset()
        shape = [dims[i][1] for i in dimids]
        vars_.append({"name": nm, "type": NC_TYPE.get(tp, tp), "typecode": tp,
                      "shape": shape, "dims": [dims[i][0] for i in dimids],
                      "units": vatts.get("units", ""), "atts": vatts,
                      "begin": begin, "vsize": vsize})
    return {"version": ver, "numrecs": numrecs, "dims": dims,
            "global": gatts, "vars": vars_}


if __name__ == "__main__":
    import sys
    for p in sys.argv[1:]:
        h = parse(p)
        print(f"\n=== {p}  (CDF-{h['version']}, numrecs={h['numrecs']}) ===")
        print("dims:", ", ".join(f"{n}={s}" for n, s in h["dims"]))
        print("global attrs:", {k: v for k, v in h["global"].items()})
        print("variables:")
        for v in h["vars"]:
            print(f"  {v['name']:16s} {v['type']:7s} {tuple(v['shape'])} "
                  f"{'x'.join(v['dims'])}  units={v['units']!r}")

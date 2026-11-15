import re, textwrap
src = open("reference/src/MGN2MECH/INCLDIR/SPC_NOCONVER.EXT").read().splitlines()

# Merge Fortran fixed-form continuation lines (leading '&').
logical = []
for ln in src:
    s = ln.rstrip()
    if s.lstrip().startswith('&'):
        logical[-1] += ' ' + s.lstrip()[1:]
    else:
        logical.append(s)

names = {}   # spca_idx(1..201) -> name
mwt   = {}   # spca_idx -> molar weight
re_spc = re.compile(r"DATA\s+SPCA_SPC\(\s*(\d+)\),\s*SPCA_MWT\(\s*\d+\)\s*/\s*'([^']*)'\s*,\s*([0-9.eE+-]+)\s*/", re.I)
for L in logical:
    m = re_spc.search(L)
    if m:
        i = int(m.group(1)); names[i] = m.group(2); mwt[i] = float(m.group(3))

# Mapping: SPCA_NAM(i), SPCA_MAP(i), MG20_NAM(i), MG20_MAP(i), EFFS_ALL(i)
mp = {}   # map entry i -> (spca_idx, mg20_idx, frac)
re_map = re.compile(
    r"SPCA_NAM\(\s*(\d+)\).*?/\s*'([^']*)'\s*,\s*(\d+)\s*,\s*'([^']*)'\s*,\s*(\d+)\s*,\s*([0-9.eE+-]+)\s*/", re.I)
for L in logical:
    m = re_map.search(L)
    if m:
        i = int(m.group(1))
        mp[i] = (int(m.group(3)), int(m.group(5)), float(m.group(6)))

print("species parsed:", len(names), " mwt:", len(mwt), " map entries:", len(mp))
assert len(names) == 201 and len(mwt) == 201 and len(mp) == 201, "count mismatch!"

# Mass-conservation sanity: speciation fractions per MEGAN class should sum to ~1
from collections import defaultdict
csum = defaultdict(float)
for i,(spca,mg20,frac) in mp.items():
    csum[mg20] += frac
CLS = ["ISOP","MBO","MT_PINE","MT_ACYC","MT_CAMP","MT_SABI","MT_AROM","NO","SQT_HR",
       "SQT_LR","MEOH","ACTO","ETOH","ACID","LVOC","OXPROD","STRESS","OTHER","CO"]
print("\n class  sum(effs_all)")
bad=0
for c in range(1,20):
    s=csum.get(c,0.0)
    flag = "" if abs(s-1.0)<1e-6 else "  <-- != 1"
    if flag: bad+=1
    print(f"  {c:2d} {CLS[c-1]:8s} {s:.5f}{flag}")
print("classes not summing to 1:", bad)

# Emit the data module
spca_names = [names[i] for i in range(1,202)]
spca_mwt   = [mwt[i]   for i in range(1,202)]
map_spca = [mp[i][0]-1 for i in range(1,202)]   # 0-based SPCA index
map_mg20 = [mp[i][1]-1 for i in range(1,202)]   # 0-based MEGAN class index
map_effs = [mp[i][2]   for i in range(1,202)]

def fmt(lst, width=6):
    body = ", ".join(repr(x) for x in lst)
    return "\n".join(textwrap.wrap(body, 96, initial_indent="    ", subsequent_indent="    "))

with open("src/megan_sc/engine/speciation_data.py","w") as f:
    f.write('"""AUTO-GENERATED from reference/src/MGN2MECH/INCLDIR/SPC_NOCONVER.EXT.\n\n')
    f.write("Do not edit by hand; regenerate with tools/parse_spca.py. Holds MEGAN's\n")
    f.write("201-species speciation table: names, molar weights (g/mol), and the map\n")
    f.write("from each of the 201 species to its parent MEGAN class (0-based) with the\n")
    f.write('speciation fraction EFFS_ALL.\n"""\n\n')
    f.write("N_SPCA = 201\n\n")
    f.write("SPCA_NAMES = [\n" + fmt([repr(n) for n in spca_names]).replace("'\"","'").replace("\"'","'") + "\n]\n")
    # simpler: write lists directly
with open("src/megan_sc/engine/speciation_data.py","w") as f:
    f.write('"""AUTO-GENERATED from reference/src/MGN2MECH/INCLDIR/SPC_NOCONVER.EXT.\n\n')
    f.write("Do not edit by hand; regenerate with tools/parse_spca.py. Holds MEGAN's\n")
    f.write("201-species speciation table: names, molar weights (g/mol), and the map from\n")
    f.write("each of the 201 speciated species to its parent MEGAN class (0-based index)\n")
    f.write('with the speciation fraction EFFS_ALL.\n"""\n\n')
    f.write("N_SPCA = 201\n\n")
    f.write("SPCA_NAMES = [\n")
    f.write(fmt(spca_names) + "\n]\n\n")
    f.write("# molar weight (g/mol) per SPCA species\n")
    f.write("SPCA_MWT = [\n" + fmt(spca_mwt) + "\n]\n\n")
    f.write("# 0-based parent MEGAN class (0..18) for each of the 201 species\n")
    f.write("MAP_MG20 = [\n" + fmt(map_mg20) + "\n]\n\n")
    f.write("# 0-based SPCA species index each map entry writes to\n")
    f.write("MAP_SPCA = [\n" + fmt(map_spca) + "\n]\n\n")
    f.write("# speciation fraction of the parent class assigned to this species\n")
    f.write("MAP_EFFS = [\n" + fmt(map_effs) + "\n]\n")
print("\nwrote src/megan_sc/engine/speciation_data.py")

# Running EIDORS v3.8 on Apple Silicon

Two obstacles, both one-off.

## 1. `eidors_var_id` has no `mexmaca64` build

EIDORS ships pre-compiled `eidors_var_id` for `mexw64`, `mexmaci64` and `mexmac`
only. On an arm64 MATLAB the startup tries to compile one and fails: MathWorks
requires a full Xcode installation for `mex` on macOS, and the Command Line
Tools alone are not accepted.

The function computes a hash used as a cache key; callers only ever compare the
returned strings. A pure-MATLAB replacement is therefore sufficient. Place this
as `eidors/arch/matlab/eidors_var_id.m`:

```matlab
function id = eidors_var_id(var)
bytes = getByteStreamFromArray(var);
d     = matlab.internal.crypto.BasicDigester('SHA256');
h     = d.computeDigest(uint8(bytes(:)));
hex   = sprintf('%02X', h);
id    = ['id_' hex(1:40)];      % keeps the original 'id_' + 40 hex format
end
```

SHA-256 truncated to 160 bits, because `BasicDigester` on R2026a accepts
`SHA256` but not `SHA1`.

Then patch `eidors/eidors_startup.m` so a missing mex falls back to the .m file
instead of aborting the whole startup:

```matlab
if ~exist(fname, 'file')
   if exist([HOMEDIR, archdir, '/eidors_var_id.m'], 'file')
      eidors_msg('STARTUP: no mex for this platform; using eidors_var_id.m', 1);
   else
      eidors_msg('STARTUP: missing a required, pre-compiled mex file: eidors_var_id', 1);
      compile_mex(HOMEDIR,archdir,ver);
   end
end
```

## 2. Netgen 6.1 is x86_64

The 6.1 build predates Apple Silicon. Install Rosetta 2:

```
softwareupdate --install-rosetta --agree-to-license
```

EIDORS invokes Netgen as a separate process, so a native arm64 MATLAB calling an
x86_64 mesher through Rosetta works. `NETGENDIR` and `PATH` must both be set
before EIDORS starts; `setup_eit_env.m` does this and then verifies it by
building a mesh.

## A trap worth knowing

`find_boundary(fmdl.elems)` returns fewer faces than the boundary EIDORS stored
when the mesh was built — 5394 against 5988 on the phantom used here — and the
missing faces render as long stray edges across the axes. Use `fmdl.boundary`
when drawing the outer surface. This is why `show_fem` produces streaked figures
on this model.

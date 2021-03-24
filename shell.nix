with (import <nixpkgs> {});

let pythonEnv = python38.withPackages (
      ps: [ps.numpy
           ps.scipy
           ps.sparse
          ]);
in mkShell rec {
  buildInputs = [pythonEnv];
}

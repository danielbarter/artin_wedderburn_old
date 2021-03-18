with (import /home/danielbarter/nixpkgs {});

let pythonEnv = python38.withPackages (
      ps: [ps.numpy
           ps.scipy
          ]);
in mkShell rec {
  buildInputs = [pythonEnv];
}

with (import <nixpkgs> {});

let pythonEnv = python38.withPackages (
      ps: [ ps.numpy
            ps.scipy
          ]);

in mkShell {

  buildInputs = [ pythonEnv
                ];

}

with (import <nixpkgs> {});

let pythonEnv = python38.withPackages (
      ps: [ ps.numpy
            ps.scipy
            ps.cupy
          ]);

in mkShell rec {

  buildInputs = [ pythonEnv
                  nvtop
                ];

  # nvtop can't find shared libraries. Should be fixed in nvtop
  shellHook = ''
    LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${pkgs.linuxPackages.nvidia_x11}/lib
   '';

}

{
  description = "A very basic flake";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/release-22.05";
    online_lib = {
        # this is a personal mirror of the EigenPositIntegration. For some reason, I cannot access the private github repo even with a token. So now I have a personal gitea instance on which this is hosted "public"
        url = "git+http://10.0.0.1:3000/aethan/EigenPositIntegration.git";
    };
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, online_lib, flake-utils }:
    let
      # wrap this in another let .. in to add the hydra job only for a single architecture
      output_set = flake-utils.lib.eachDefaultSystem (system:
        let
            pkgs = nixpkgs.legacyPackages.${system};
            python = pkgs.python3.withPackages(p: with p; [numpy]);
        in
        rec {
            packages = flake-utils.lib.flattenTree {
                cim = pkgs.gcc10Stdenv.mkDerivation {
                    name = "CIM";
                    src = ./.;

                    nativeBuildInputs = [pkgs.cmake];

                    buildInputs = [
                        online_lib.packages.${system}.eigen
                        online_lib.packages.${system}.universal
                        online_lib.packages.${system}.eigen_universal_integration
                        pkgs.llvmPackages.openmp
                    ];

                    checkPhase = ''
                        ctest
                    '';

                    cmakeFlags = [
                        "-DTESTS=ON"
                        "-DDEBUG=OFF"
                    ];

                    doCheck = true;
                };
            };

            defaultPackage = packages.cim;

            devShell = pkgs.mkShell {
                buildInputs = [
                    online_lib.packages.${system}.universal
                    online_lib.packages.${system}.eigen
                    pkgs.cmake
                    python
                    # pkgs.gdbgui
                ];

                shellHook = ''
                    function run_cmake_build() {
                        cd project
                        cmake -B /mnt/RamDisk/build -DDEBUG=ON -DTESTS=ON -DSOLVERS=ON -DCIM=ON
                        cmake --build /mnt/RamDisk/build -j1
                        cd ..
                    }

                    function build_and_run() {
                        run_cmake_build
                        /mnt/RamDisk/build/ShermanLoad
                    }
                    function build_and_test() {
                        run_cmake_build
                        pushd /mnt/RamDisk/build
                        ctest
                        popd
                    }
                '';
            };


        }
    );
    in
        output_set // { hydraJobs.build."aarch64-linux" = output_set.defaultPackage."aarch64-linux"; };
    }

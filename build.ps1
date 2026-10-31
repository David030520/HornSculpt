param(
    [Parameter(Mandatory = $true)][string]$QtRoot,
    [Parameter(Mandatory = $true)][string]$CgalRoot,
    [Parameter(Mandatory = $true)][string]$BoostRoot,
    [Parameter(Mandatory = $true)][string]$OcctRoot,
    [Parameter(Mandatory = $true)][string]$OcctRuntimeDirectory,
    [string]$CudaRoot = $env:CUDA_PATH,
    [string]$BuildDirectory = 'build',
    [ValidateSet('Release', 'Debug', 'RelWithDebInfo')][string]$Configuration = 'Release',
    [switch]$Test
)

$ErrorActionPreference = 'Stop'
$sourceDirectory = $PSScriptRoot
if (-not [System.IO.Path]::IsPathRooted($BuildDirectory)) {
    $BuildDirectory = Join-Path $sourceDirectory $BuildDirectory
}
$qtDirectory = (Resolve-Path -LiteralPath $QtRoot).Path
$cgalDirectory = (Resolve-Path -LiteralPath $CgalRoot).Path
$boostDirectory = (Resolve-Path -LiteralPath $BoostRoot).Path
$cudaDirectory = (Resolve-Path -LiteralPath $CudaRoot).Path
$occtDirectory = (Resolve-Path -LiteralPath $OcctRoot).Path
$occtRuntime = (Resolve-Path -LiteralPath $OcctRuntimeDirectory).Path
& cmake -S $sourceDirectory -B $BuildDirectory -G 'Visual Studio 17 2022' -A x64 `
    -T "cuda=$cudaDirectory" "-DCMAKE_PREFIX_PATH=$qtDirectory;$cgalDirectory;$boostDirectory" `
    -DCGAL_DISABLE_GMP=ON -DCMAKE_CUDA_ARCHITECTURES=native -DBUILD_TESTING=ON `
    "-DHORNSCULPT_OCCT_ROOT=$occtDirectory" "-DHORNSCULPT_OCCT_RUNTIME_DIR=$occtRuntime"
if ($LASTEXITCODE -ne 0) { throw 'CMake-Konfiguration fehlgeschlagen.' }
& cmake --build $BuildDirectory --config $Configuration --parallel
if ($LASTEXITCODE -ne 0) { throw 'C++-Build fehlgeschlagen.' }
if ($Test) {
    $env:PATH = (Join-Path $cudaDirectory 'bin/x64') + ';' + $env:PATH
    & ctest --test-dir $BuildDirectory -C $Configuration --output-on-failure
    if ($LASTEXITCODE -ne 0) { throw 'Tests fehlgeschlagen.' }
}

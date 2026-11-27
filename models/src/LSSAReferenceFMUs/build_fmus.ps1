# Builds the Windows binaries (x86_64-windows) of the LS-SA NBouncingBalls FMUs for all ball counts and adds them to
#   models/bin/LSSAReferenceFMUs/<Version>/3.0/NBouncingBalls_<N>.fmu
# Binaries of other platforms already in the FMU (x86_64-linux from build_fmus.sh) are kept. If an existing FMU has a
# different modelDescription.xml, its binaries were built from different sources: the script stops.
#
# usage (from any directory):
#   powershell -ExecutionPolicy Bypass -File build_fmus.ps1 [-Version 0.2.0] [-Counts 1,3,10,30,100,300,1000]
#
# Uses `cmake` from PATH, or the CMake bundled with the latest Visual Studio installation.
# All FMUs of one version must be built with the same compiler / flags (they are used for benchmarks).
#
# Variants (e.g. for benchmarking without the plan cache) go to a separate build and output directory:
#   build_fmus.ps1 -BuildDir build-nocache -OutDir C:\tmp\nocache -CMakeArgs "-DLSSA_PLAN_CACHE=OFF"
#
# Coupled variant (NCoupledBouncingBalls_<N>_K<K>.fmu): -Neighbors lists the neighbors per side K; K is capped at
# N - 1 (-1: dense, K = N - 1), N = 1 (no neighbors) and duplicates are skipped:
#   build_fmus.ps1 -Counts 3,10,30,100,300,1000 -Neighbors 1,8,-1
param(
    [string]$Version = "0.2.0",
    [int[]]$Counts = @(1, 3, 10, 30, 100, 300, 1000),
    [string]$Generator = "Visual Studio 17 2022",
    [string]$CMake = "",
    [string]$BuildDir = "build",
    [string]$OutDir = "",
    [string[]]$CMakeArgs = @(),
    [int[]]$Neighbors = @()
)

$ErrorActionPreference = "Stop"

$src = $PSScriptRoot
$out = if ($OutDir) { $OutDir } else { Join-Path $src "..\..\bin\LSSAReferenceFMUs\$Version\3.0" }
$out = [System.IO.Path]::GetFullPath($out)

if (-not $CMake) {
    $cmd = Get-Command cmake -ErrorAction SilentlyContinue
    if ($cmd) { $CMake = $cmd.Source }
}
if (-not $CMake) {
    $vswhere = Join-Path ${env:ProgramFiles(x86)} "Microsoft Visual Studio\Installer\vswhere.exe"
    $vs = & $vswhere -latest -property installationPath
    $CMake = Join-Path $vs "Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe"
}
if (-not (Test-Path $CMake)) { throw "CMake not found." }

$platformArgs = @()
if ($Generator -like "Visual Studio*") { $platformArgs = @("-A", "x64") }

New-Item -ItemType Directory -Force $out | Out-Null

Add-Type -AssemblyName System.IO.Compression.FileSystem
$platform = "x86_64-windows"

# adds the content of the staging directory (modelDescription.xml, extra/, binaries/<platform>/) to the FMU $target
function Merge-FmuContent([string]$staging, [string]$target) {
    $tmp = Join-Path ([System.IO.Path]::GetTempPath()) ([System.IO.Path]::GetRandomFileName())
    New-Item -ItemType Directory $tmp | Out-Null
    try {
        if (Test-Path $target) {
            [System.IO.Compression.ZipFile]::ExtractToDirectory($target, $tmp)
            $old = (Get-Content -Raw (Join-Path $tmp "modelDescription.xml")) -replace "`r", ""
            $new = (Get-Content -Raw (Join-Path $staging "modelDescription.xml")) -replace "`r", ""
            if ($old -ne $new) {
                throw "$target has a different modelDescription.xml (built from different sources); rebuild all platforms or remove the FMU first"
            }
            $old = Join-Path $tmp "binaries\$platform"
            if (Test-Path $old) { Remove-Item -Recurse -Force $old }
        }
        Copy-Item -Recurse -Force (Join-Path $staging "*") $tmp
        $entries = Get-ChildItem -Force $tmp | ForEach-Object { $_.Name }
        Push-Location $tmp
        try {
            & $CMake -E tar cf "$target.tmp" --format=zip -- @entries
            if ($LASTEXITCODE -ne 0) { throw "packing $target failed" }
        } finally { Pop-Location }
        Move-Item -Force "$target.tmp" $target
    } finally {
        Remove-Item -Recurse -Force $tmp
    }
}

# (N, K) pairs to build: K = 0 is the uncoupled model
$builds = @()
foreach ($n in $Counts) {
    if ($Neighbors.Count -eq 0) {
        $builds += , @($n, 0)
        continue
    }
    foreach ($k in $Neighbors) {
        $keff = if ($k -lt 0 -or $k -gt $n - 1) { $n - 1 } else { $k }
        if ($keff -lt 1) { continue }
        if (-not ($builds | Where-Object { $_[0] -eq $n -and $_[1] -eq $keff })) { $builds += , @($n, $keff) }
    }
}

foreach ($nk in $builds) {
    $n = $nk[0]
    $k = $nk[1]
    $name = if ($k -eq 0) { "NBouncingBalls_$n" } else { "NCoupledBouncingBalls_${n}_K$k" }
    $build = if ($k -eq 0) { Join-Path $src "$BuildDir\n$n" } else { Join-Path $src "$BuildDir\n${n}_k$k" }
    Write-Host "=== $name ==="

    & $CMake -S $src -B $build -G $Generator @platformArgs "-DNUM_BALLS=$n" "-DNUM_NEIGHBORS=$k" @CMakeArgs
    if ($LASTEXITCODE -ne 0) { throw "CMake configure failed for $name" }

    & $CMake --build $build --config Release
    if ($LASTEXITCODE -ne 0) { throw "Build failed for $name" }

    $target = Join-Path $out "$name.fmu"
    Merge-FmuContent (Join-Path $build "fmu_staging") $target
    Write-Host "-> $target"
}

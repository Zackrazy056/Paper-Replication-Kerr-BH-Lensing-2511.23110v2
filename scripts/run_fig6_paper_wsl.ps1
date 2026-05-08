param(
    [string]$Python = "",
    [string]$Distro = "Ubuntu",
    [string]$PythonVersion = "312",
    [string]$Abi = "cp312",
    [string]$Out = "results/Fig6_paper_standard.png",
    [string]$Summary = "results/Fig6_paper_standard_summary.json",
    [int]$NFreq = 2400,
    [switch]$SkipDownload
)

$ErrorActionPreference = "Stop"

$Root = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
Set-Location $Root

function Get-WslPath {
    param([string]$WindowsPath)

    $output = & wsl -d $Distro -- wslpath -a "$WindowsPath"
    $path = ($output | Where-Object { $_ -like "/*" } | Select-Object -Last 1)
    if ($path) {
        return $path.Trim()
    }
    return ""
}

if (-not $Python) {
    $parentPython = Resolve-Path (Join-Path $Root "..\.venv\Scripts\python.exe") -ErrorAction SilentlyContinue
    $localPython = Resolve-Path (Join-Path $Root ".venv\Scripts\python.exe") -ErrorAction SilentlyContinue
    if ($parentPython) {
        $Python = $parentPython.Path
    } elseif ($localPython) {
        $Python = $localPython.Path
    } else {
        $Python = "python"
    }
}

if (-not $SkipDownload) {
    & $Python -m pip download `
        --dest wheelhouse-paper-wsl `
        --index-url https://pypi.org/simple `
        --platform manylinux_2_28_x86_64 `
        --platform manylinux2014_x86_64 `
        --platform manylinux_2_17_x86_64 `
        --python-version $PythonVersion `
        --implementation cp `
        --abi $Abi `
        --only-binary=:all: `
        -r requirements-paper.txt
    if ($LASTEXITCODE -ne 0) {
        throw "Linux wheel download failed."
    }
}

$wslRoot = Get-WslPath $Root
if ($LASTEXITCODE -ne 0 -or -not $wslRoot) {
    throw "Could not translate project path for WSL distro '$Distro'."
}

$wslRootQuoted = "'" + $wslRoot.Replace("'", "'\''") + "'"
$outQuoted = "'" + $Out.Replace("'", "'\''") + "'"
$summaryQuoted = "'" + $Summary.Replace("'", "'\''") + "'"

$bash = @"
set -euo pipefail
cd $wslRootQuoted

if [ ! -x .venv-paper-wsl/bin/python ]; then
  rm -rf .venv-paper-wsl
  if python3 -m venv .venv-paper-wsl; then
    :
  else
    rm -rf .venv-paper-wsl
    if python3 -m virtualenv .venv-paper-wsl; then
      :
    else
      echo "Could not create .venv-paper-wsl." >&2
      echo "Install python3-venv in WSL, or install virtualenv, then rerun this script." >&2
      exit 2
    fi
  fi
fi

. .venv-paper-wsl/bin/activate
python -m pip install --quiet --no-index --find-links wheelhouse-paper-wsl -r requirements-paper.txt
python main_fig6_paper_standard.py --check-deps
MPLBACKEND=Agg python main_fig6_paper_standard.py --strict-sr --n-freq $NFreq --out $outQuoted --summary $summaryQuoted
"@

$tmpScript = Join-Path $Root ".tmp_fig6_paper_wsl.sh"
$bashLf = $bash.Replace("`r`n", "`n")
[System.IO.File]::WriteAllText($tmpScript, $bashLf, [System.Text.Encoding]::ASCII)

$wslScript = Get-WslPath $tmpScript
if ($LASTEXITCODE -ne 0 -or -not $wslScript) {
    Remove-Item -LiteralPath $tmpScript -ErrorAction SilentlyContinue
    throw "Could not translate temporary script path for WSL distro '$Distro'."
}

& wsl -d $Distro -- bash "$wslScript"
$exitCode = $LASTEXITCODE
Remove-Item -LiteralPath $tmpScript -ErrorAction SilentlyContinue

if ($exitCode -ne 0) {
    throw "WSL paper-standard Fig.6 run failed."
}

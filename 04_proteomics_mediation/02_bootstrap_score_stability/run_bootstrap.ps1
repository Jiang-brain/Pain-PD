param(
    [string]$Rscript = 'Rscript',
    [ValidateRange(1, 16)][int]$Workers = 4
)
$ErrorActionPreference = 'Stop'
$codeRoot = (Resolve-Path (Join-Path $PSScriptRoot '../..')).Path
$rCommand = (Get-Command $Rscript -ErrorAction Stop).Source
$resultsRoot = if ($env:PD_PAIN_RESULTS_DIR) { $env:PD_PAIN_RESULTS_DIR } else { Join-Path $codeRoot 'results' }
$output = Join-Path $resultsRoot 'sensitivity_nonproteomic_weights\bootstrap1000'
New-Item -ItemType Directory -Path $output -Force | Out-Null
Push-Location $codeRoot
try {
    & $rCommand '04_proteomics_mediation/02_bootstrap_score_stability/01_prepare_bootstrap_inputs.R'
    if ($LASTEXITCODE -ne 0) { throw 'Bootstrap input preparation failed.' }
    $jobs = foreach ($shard in 1..$Workers) {
        Start-Process -FilePath $rCommand -WorkingDirectory $codeRoot -WindowStyle Hidden -PassThru `
            -ArgumentList @('04_proteomics_mediation/02_bootstrap_score_stability/02_run_bootstrap.R', "$shard", "$Workers") `
            -RedirectStandardOutput (Join-Path $output "shard_${shard}_stdout.log") `
            -RedirectStandardError (Join-Path $output "shard_${shard}_stderr.log")
    }
    $jobs | Select-Object Id, ProcessName
    $jobs | Wait-Process
    $exitCodes = $jobs | ForEach-Object { $_.Refresh(); $_.ExitCode }
    if ($exitCodes | Where-Object { $_ -ne 0 }) { throw "One or more bootstrap shards failed: $exitCodes" }
    foreach ($step in @(
        @('04_proteomics_mediation/02_bootstrap_score_stability/02_run_bootstrap.R', '--collect-only'),
        @('04_proteomics_mediation/02_bootstrap_score_stability/03_score_correlations.R'),
        @('04_proteomics_mediation/02_bootstrap_score_stability/04_validate_and_report.R')
    )) {
        & $rCommand @step
        if ($LASTEXITCODE -ne 0) { throw "Bootstrap step failed: $step" }
    }
} finally {
    Pop-Location
}

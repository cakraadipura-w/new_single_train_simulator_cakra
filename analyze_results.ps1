$df = Import-Csv "experiment_results/IS02/E5_results.csv"
$metrics = @("hv", "igd", "feasible", "E_best_kWh", "energy_saving_pct", "runtime_s")
$grouped = $df | Group-Object config
$summary = @()
foreach ($group in $grouped) {
    $item = [PSCustomObject]@{ config = $group.Name }
    foreach ($m in $metrics) {
        $values = $group.Group | ForEach-Object { [double]$_.$m }
        $avg = ($values | Measure-Object -Average).Average
        $sumSq = 0
        foreach ($v in $values) { $sumSq += [Math]::Pow($v - $avg, 2) }
        $std = [Math]::Sqrt($sumSq / ($values.Count - 1))
        $item | Add-Member -MemberType NoteProperty -Name "$m`_mean" -Value $avg
        $item | Add-Member -MemberType NoteProperty -Name "$m`_std" -Value $std
    }
    $summary += $item
}
$summary | Format-Table -AutoSize
function Get-Mean { param($conf) $summary | Where-Object { $_.config -eq $conf } }
$c = Get-Mean "C"
$c_llm = Get-Mean "C_LLM"
$d = Get-Mean "D"
$d_llm = Get-Mean "D_LLM"
Write-Host "`nDeltas (C_LLM - C):"
foreach ($m in $metrics) {
    if ($c_llm -and $c) {
        $delta = $c_llm."$($m)_mean" - $c."$($m)_mean"
        Write-Host "${m}: $delta"
    }
}
Write-Host "`nDeltas (D_LLM - D):"
foreach ($m in $metrics) {
    if ($d_llm -and $d) {
        $delta = $d_llm."$($m)_mean" - $d."$($m)_mean"
        Write-Host "${m}: $delta"
    }
}

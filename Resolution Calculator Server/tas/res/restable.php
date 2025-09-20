<?php
rep($abc,3);
rep($alphabetagamma,3);
rep($hcoll,4);
rep($vcoll,4);

$content = <<<EOT
#a b c alpha beta gamma
$abc $alphabetagamma
#horizontal collimation 1 2 3 4
$hcoll
#vertical collimation 1 2 3 4
$vcoll
#monochromator mosaic
$monochromator_mosaic
#analyzer mosaic
$analyzer_mosaic
#sample horizontal mosaic
$sampleh_mosaic
#sample vertical mosaic
$samplev_mosaic
#fixed energy
$fixed_energy
#use fixed energy
$fixed_energy_choice
#monochromator tau
$monochromator_tau
#analyzer tau
$analyzer_tau
#orientation vector 1
$orientation_1
#orientation vector 2
$orientation_2
#h k l
$h
$k
$l
#omega (energy transfer)
$omega
EOT;

$infile = tempnam('/tmp','RIN');
$fid = fopen($infile, 'w');
fputs($fid,$content);
fclose($fid);

$src = getcwd();
$cmd = "$PYTHON $src/rescalc.py --table $infile 2>&1";
chdir("/tmp");
putenv("MPLCONFIGDIR=/tmp");
system($cmd, $status);
if ($status != 0) {
   echo("input file:<pre>\n");
   readfile($infile);
   echo("</pre>\n");
}
unlink($infile);
?>

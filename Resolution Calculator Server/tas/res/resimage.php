<?php
include 'resdef.php';
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

#system("rm /tmp/RIN*");
$infile = tempnam('/tmp','RIN');
$fid = fopen($infile, 'w');
fputs($fid,$content);
fclose($fid);

$src = getcwd();
#$cmd = "$PYTHON $src/rescalc.py --plot $infile 2>&1 >/dev/null";
$cmd = "$PYTHON $src/rescalc.py --plot $infile 2>&1 >/tmp/rescalc.log";
chdir("/tmp");
putenv("MPLCONFIGDIR=/tmp");
error_log('ERROR TEST');
#system("$PYTHON -V");
system($cmd, $status);
if ($status == 0) {
  $path = '/tmp/resout.png';
  $img = getimagesize($path);
  header('Content-Type: '.$img['mime']);
  header('Content-Disposition: inline; filename="resolution.png";');
  header('Content-Length: '.filesize($path));
  readfile($path);
  #unlink($path);
} else {
  echo("<html><body>\n");
  echo("$cmd<p>\n");
  echo("Could not compute resolution from:<pre>\n");
  readfile($infile);
  echo("</pre></body></html>");
}
#unlink($infile);
?>

<?php
$taus = array('pg(002)','pg(004)','ge(111)','ge(220)','ge(311)','be(002)','pg(110)','2axis');
$abc = '6.28';
$alphabetagamma = '90';
$hcoll = '40';
$vcoll = '120';
$monochromator_mosaic = '30';
$monochromator_tau = $taus[0];
$analyzer_mosaic = '30';
$analyzer_tau = $taus[0];
$sampleh_mosaic = '30';
$samplev_mosaic = '30';
$orientation_1 = '1 0 0';
$orientation_2 = '0 1 0';
$sampleh = '30';
$samplev = '30';
$fixed_energy = '14.7';
$use_fixed_energy = '1';
$h = '1 1';
$k = '0 0';
$l = '0 0';
$omega = '0 0';


# Build a select option list with one value selected
function select_box($target, $options, $size=1) {
  $value = $GLOBALS[$target];
  echo ("<select name='$target' id='$target' size='$size'>");
  foreach ($options as $opt) {
    if ($value == $opt) {
     echo ("<option selected>$opt</option>");
    } else {
     echo ("<option>$opt</option>");
    }
  }
  echo ("</select>");
}
// Build a numeric input box
function float_box($target, $size=1) {
  $value = $GLOBALS[$target];
  echo ("<input type='text' name='$target' id='$target' value='$value' />");
}
// Build a check box
function check_box($target) {
  $value = $GLOBALS[$target];
  if ($value == 1) {
    $checked = "checked";
  } else {
    $checked = "";
  }
  echo ("<input type='checkbox' name='$target' id='$target' value='1' $checked />");
}
function label($target, $label) {
  echo ("<label for='$target'>$label</label>");
}
function join_args() {
  $req='?';
  foreach (func_get_args() as $arg) {
    $req = $req.$arg.'='.$GLOBALS[$arg].'&';
  }
  return substr($req, 0, -1);
}
?>

<html>
<head> <title>Triple-Axis resolution calculator</title> </head>
<body>

<form action='rescalc.php' method='get'>
<table><tr>
<td><label for='abc'>a,b,c</label></td> 
<td><?php float_box('abc', 3); ?></td> </tr><tr>
<td><label for='alphabetagamma'>&alpha;,&beta;,&gamma;</label></td> 
<td><?php float_box('alphabetagamma', 3); ?></td>
</tr></table>

<hr />

<table><tr>
<td><label for='hcoll'>Horizontal collimator 1,2,3,4</td>
<td><?php float_box('hcoll', 4); ?></td><td>mm</td>
</tr><tr>
<td><label for='vcoll'>Vertical collimator 1,2,3,4</td>
<td><?php float_box('vcoll', 4); ?></td><td>mm</td>
</tr></table>

<hr />

<table><tr>
<td><label for='monochromator_mosaic'>Monochromator mosaic</td>
<td><?php float_box('monochromator_mosaic'); ?></td>
<td><label for='monochromator_tau'>&tau;</td>
<td><?php select_box('monochromator_tau', $taus); ?></td>
</tr><tr>
<td><label for='analyzer_mosaic'>Analyzer mosaic</td>
<td><?php float_box('analyzer_mosaic'); ?></td>
<td><label for='analyzer_tau'>&tau;</td>
<td><?php select_box('analyzer_tau', $taus); ?></td>
</tr><tr>
<td><label for='sampleh_mosaic'>Sample horz mosaic</td>
<td><?php float_box('sampleh_mosaic'); ?></td>
</tr><tr>
<td><label for='samplev_mosaic'>Sample vert mosaic</td>
<td><?php float_box('samplev_mosaic'); ?></td>
</tr></table>
<?php check_box('use_fixed_energy'); ?> 
<label for='fixed_energy'>Fixed energy</label> 
<?php float_box('fixed_energy'); ?>  <br />

<hr />

<label for='orientation_1'>Orientation 1</label>
<?php float_box('orientation_1', 3); ?>
<label for='orientation_2'>2</label>
<?php float_box('orientation_2', 3); ?>
<br />

<hr />
h <?php float_box('h', 2); ?>
k <?php float_box('k', 2); ?>
l <?php float_box('l', 2); ?>
<br />
&omega; <?php float_box('omega', 2); ?> (energy transfer)
<br />
<input type="submit" name="action" value="Plot">
</form>

<hr />
<?php
$request=join_args('abc','alphabetagamma','hcoll','vcoll',
      'monochromator_mosaic','monochromator_tau',
      'analyzer_mosaic','analyzer_tau',
      'sampleh_mosaic','samplev_mosaic',
      'orientation_1','orientation_2',
      'fixed_energy','use_fixed_energy',
      'h','k','l','omega');
echo("<img src='resimage.php$request' alt="resolution plot" />");
?>

</body>
</html>

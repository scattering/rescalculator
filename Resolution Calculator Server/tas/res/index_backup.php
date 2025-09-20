<?php
include 'resdef.php';
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
<head> 
<title>Triple-Axis resolution calculator</title>

<script src="jquery-1.1.3.1.pack.js" type="text/javascript"></script>
<script src="jquery.tabs.pack.js" type="text/javascript"></script>
<script type="text/javascript" src="jquery.form.js"></script> 

<link rel="stylesheet" href="site.css" type="text/css" media="print, projection, screen">
<link rel="stylesheet" href="jquery.tabs.css" type="text/css" media="print, projection, screen">
<!-- Additional IE/Win specific style sheet (Conditional Comments) -->
<!--[if lte IE 7]>
<link rel="stylesheet" href="jquery.tabs-ie.css" type="text/css" media="projection, screen">
<![endif]-->


<script type="text/javascript">

$(function() {
  $('#input_menu_list').show()
  $('#output_menu_list').show()
  $('#graph').append("<p id='loading' style='display: none;'></p>");
  $('#input_menu').tabs({ fxAutoHeight: true });
  $('#output_menu').tabs();
  //$('#input_form').ajaxForm();
  $('#plot').load(function(){
    $('#loading').hide();
    $('#plot').show();
  });
  $('#plot').error(function(){
    $('#loading').text('Could not plot resolution.  See table for details');
  });
  $('#plot_button').click(function (){
    var q=$('#input_form').formSerialize(),
        plot=$('#plot'), 
        table=$('#table'), 
        loading=$('#loading');
    plot.hide();
    loading.text('Calculating resolution...');
    loading.show();
    var src = plot.attr("src");
    plot.attr("src","resimage.php?"+q);
    if (plot.attr("src") == src) { 
        loading.hide();
        plot.show();
    } else {
        table.html('<p>Calculating resolution...</p>');
        $.post("restableajax.php",q,function(data){
              table.html(data);
            }, "html");
    }
    return false;
  });
  //$(window).unload(function() { alert("unloading"); return false; } );
});
</script>
</head>

<body>

<form id="input_form" action='index.php' method='post'>

<div id="input_menu">
<ul class="anchors" id="input_menu_list" style="display: none;">
<li><a href="#collimations">Collimations</a></li>
<li><a href="#hklw">H-K-L-&omega;</a></li>
<li><a href="#energy">Energy</a></li>
<li><a href="#mosaic">Mosaic</a></li>
<li><a href="#lattice">Lattice</a></li>
<li><a href="#orientation">Orientation</a></li>
</ul>

<div id='orientation' class='anchor'>
<table><tr>
<td><label for='orientation_1'>Orientation 1</label></td>
<td><?php float_box('orientation_1', 3); ?></td>
</tr><tr>
<td><label for='orientation_2'>Orientation 2</label></td>
<td><?php float_box('orientation_2', 3); ?></td>
</tr></table>
</div>

<div id='lattice' class='anchor'>
<table><tr>
<td><label for='abc'>a b c</label></td> 
<td><?php float_box('abc', 3); ?></td> 
<td style="line-height: 1ex;">&#8491; (e.g., <em>6.28</em> or <em>3 4 5</em>)</td> <!-- Firefox on windows not rendering A properly -->
</tr><tr>
<td><label for='alphabetagamma'>&alpha; &beta; &gamma;</label></td> 
<td><?php float_box('alphabetagamma', 3); ?></td>
<td>degrees (e.g., <em>90</em> or <em>90 90 120</em>)</td>
</tr></table>
</div>

<div id='mosaic' class='anchor'>
<table><tr>
<td><label for='monochromator_mosaic'>Monochromator mosaic</td>
<td><?php float_box('monochromator_mosaic'); ?></td>
<td>min</td>
<td><label for='monochromator_tau'>&tau;</td>
<td><?php select_box('monochromator_tau', $mono_taus); ?></td>
</tr><tr>
<td><label for='analyzer_mosaic'>Analyzer mosaic</td>
<td><?php float_box('analyzer_mosaic'); ?></td>
<td>min</td>
<td><label for='analyzer_tau'>&tau;</td>
<td><?php select_box('analyzer_tau', $anal_taus); ?></td>
</tr><tr>
<td><label for='sampleh_mosaic'>Sample horz mosaic</td>
<td><?php float_box('sampleh_mosaic'); ?></td>
<td>min</td>
</tr><tr>
<td><label for='samplev_mosaic'>Sample vert mosaic</td>
<td><?php float_box('samplev_mosaic'); ?></td>
<td>min</td>
</tr></table>
</div>

<div id='energy' class='anchor'>
<table><tr>
<td><label for='fixed_energy'>Fixed energy</label> </td>
<td><?php select_box('fixed_energy_choice',array('Ei','Ef')) ?></td>
<td><?php float_box('fixed_energy'); ?></td>
<td>meV</td>
</tr></table>
</div>

<div id='collimations' class='anchor'>
<table><tr>
<td><label for='hcoll'>Horizontal collimator</td>
<td><?php float_box('hcoll', 4); ?></td>
<td>min (e.g., <em>40</em> or <em>40 47 40 200</em>)</td>
</tr><tr>
<td><label for='vcoll'>Vertical collimator</td>
<td><?php float_box('vcoll', 4); ?></td>
<td>min (e.g., <em>120</em> or <em>120 120 120 120</em>)</td>
</tr></table>
</div>


<div id='hklw' class='anchor'>
Enter one value of h, k, l, &omega; for each point:
<table><tr>
<td>h</td> <td><?php float_box('h', 0); ?></td><td>(e.g., <em>1</em> or <em>1 2 3 4 &#0133;</em>)
</tr><tr>
<td>k</td> <td><?php float_box('k', 0); ?></td>
</tr><tr>
<td>l</td> <td><?php float_box('l', 0); ?></td>
</tr><tr>
<td>&omega;</td> <td> <?php float_box('omega', 2); ?> </td>
<td>meV</td><td>(energy transfer)</td>
</tr><tr>
</tr></table>
</div>
</div>

<input id="plot_button" type="submit" name="action" value="Calculate">
</form>

<div id="output_menu">
<ul class="anchors" id="output_menu_list" style="display: none">
<li><a href="#graph">Graph</a></li>
<li><a href="#table">Table</a></li>
</ul>

<div id="graph" class="anchor">
<?php
if (isset($_REQUEST['action'])) {
  $request=join_args('abc','alphabetagamma','hcoll','vcoll',
      'monochromator_mosaic','monochromator_tau',
      'analyzer_mosaic','analyzer_tau',
      'sampleh_mosaic','samplev_mosaic',
      'orientation_1','orientation_2',
      'fixed_energy','use_fixed_energy',
      'h','k','l','omega');
  echo("<img id='plot' src='resimage.php$request' alt='resolution plot' />");
} else {
  echo("<img id='plot' src='' alt='resolution plot' style='display: none;' />");
}
?>
</div>

<div id="table" class="anchor">
<?php 
if (isset($_REQUEST['action'])) { include 'restable.php'; }
?>
</div>
</div>

</body>
</html>

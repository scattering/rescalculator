<?php

#$PYTHON="/Library/Frameworks/Python.framework/Versions/Current/bin/python";
#$PYTHON="/usr/bin/python";
$PYTHON="/home/ylem/.conda/envs/drneutron/bin/python";

# assign the value from the request or set a default.
function def($name,$value) {
   if (isset($_REQUEST[$name])) { 
      $value = $_REQUEST[$name]; 
   }
   $GLOBALS[$name] = $value;
}
# repeat the value if only one value is given.
function rep(&$name, $count) {
   if (strpos($name, ' ')) {
   } else {
      $base = $name;
      for ($i = 1; $i < $count; $i++) {
          $name = $name.' '.$base;
      }
   }
}

$mono_taus = array('pg(002)','pg(004)','ge(111)','ge(220)',
              'ge(311)','be(002)','pg(110)');
$anal_taus = array('pg(002)','pg(004)','ge(111)','ge(220)',
              'ge(311)','be(002)','pg(110)','2axis');
def('abc', '6.28 6.28 6.28');
def('alphabetagamma', '90');
def('hcoll', '40');
def('vcoll', '120');
def('monochromator_mosaic', '30');
def('monochromator_tau', $mono_taus[0]);
def('analyzer_mosaic', '30');
def('analyzer_tau', $anal_taus[0]);
def('sampleh_mosaic', '30');
def('samplev_mosaic', '30');
def('orientation_1', '1 0 0');
def('orientation_2', '0 1 0');
def('fixed_energy', '14.7');
def('fixed_energy_choice', 'Ef');
def('h', '1');
def('k', '0');
def('l', '0');
def('omega', '0');

?>

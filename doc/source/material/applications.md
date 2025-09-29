# Applications that Use SLEPc

The list below shows publications on applications that use SLEPc. If you want your application to be included in this list please [contact us](../contact/index).  This list helps us demonstrate to our sponsors the usefulness of our work, so additional entries help us to continue developing, improving, and supporting SLEPc.

<script src="../_static/js/apps.js"></script>

<div id="sortedbyapp">
<p>
Applications are currently grouped by <b>application area</b>.
<a href="#" id="sortedbyprob-show" onclick="toggleview('sortedbyapp','sortedbyprob');return false;">Click here if you want them grouped by problem type.</a>
</p>
<script>
var appareas = [];
for (var i=0;i<14;i++) {
  appareas[i] = { index: i, count: countapps(i) };
}
appareas.sort(function (a,b) {
  return b.count-a.count;
});
document.write('<table class=apps cellpadding=0px cellspacing=0px>');
document.write('<tr><th>Area</th><th>Count</th></tr>');
for (var j=0;j<14;j++) {
  var i = appareas[j].index;
  document.write('<tr><td><a href="#'+categword(i)+'">'+categ[i]+'</a></td>');
  document.write('<td>'+countapps(i)+'</td></tr>');
}
document.write('<tr><td><i>Total</i></td>');
document.write('<td><b>'+apps.length+'</b></td></tr>');
document.write('</table>');
for (var j=0;j<14;j++) {
  var i = appareas[j].index;
  document.write('<h2 id="'+categword(i)+'">'+categ[i]+'</h2>');
  document.write('<ol>');
  listapps(i);
  document.write('</ol>');
}
</script>
</div>

<div id="sortedbyprob" style="display: none">
<p>
Applications are currently grouped by <b>problem type</b>.
<a href="#" id="sortedbyapp-show" onclick="toggleview('sortedbyprob','sortedbyapp');return false;">Click here if you want them grouped by application area.</a>
</p>
<script>
var probcount = [];
for (var i=0;i<7;i++) {
  probcount[i] = { index: i, count: countprob(i) };
}
probcount.sort(function (a,b) {
  return b.count-a.count;
});
document.write('<table class=apps cellpadding=0px cellspacing=0px>');
document.write('<tr><th>Problem type</th><th>Count</th></tr>');
for (var j=0;j<7;j++) {
  var i = probcount[j].index;
  document.write('<tr><td><a href="#'+probword(i)+'">'+probtype[i]+'</a></td>');
  document.write('<td>'+countprob(i)+'</td></tr>');
}
document.write('<tr><td><i>Total</i></td>');
document.write('<td><b>'+apps.length+'</b></td></tr>');
document.write('</table>');
for (var j=0;j<7;j++) {
  var i = probcount[j].index;
  document.write('<h2 id="'+probword(i)+'">'+probtype[i]+'</h2>');
  document.write('<ol>');
  listappsbyprob(i);
  document.write('</ol>');
}
</script>
</div>


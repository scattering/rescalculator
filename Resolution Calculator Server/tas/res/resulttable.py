import numpy

def render(table):

    names = table[::2]
    vals = table[1::2] 
    print('<table class="data">')
    print('<thead><tr class="odd">', end=' ')
    for n in names:
        print('<th>',n,'</th>', end=' ')
    print('</tr></thead>')

    for i,row in enumerate(zip(*vals)):
        if i%2: htmlclass='class="odd"' 
        else: htmlclass=''
        print("<tr %s>"%htmlclass)
        for col in row:
            if numpy.isscalar(col):
                print("<td>%.4f</td>"%col, end=' ')
            elif isinstance(col, str):
                print("<td>%s</td>"%col, end=' ')
            elif isinstance(col,numpy.ndarray):
                print('<td class="embed">')
                render_array(col)
                print("</td>")
            else:
                print("<td></td>", end=' ')
        print("</tr>") 

    print("</table>")

popup_num = 0
def render_array(A):
    global popup_num
    popup_num += 1
    opt = dict(id="popup%d"%popup_num)
#    print '''<style type="text/css">
##%(id)s a:hover .popup { visibility: visible; }
##%(id)s a .popup { visibility: hidden; }
#</style>'''%opt
#    print '<a href="#popup%d" id="popup%d">[]<span class="popup">'%(popup_num,popup_num)
    print('<table class="data,embed">', end=' ')
    for row in A:
        print("<tr>", end=' ')
        for col in row:
            print("<td>%.4f</td>"%col, end=' ')
        print("</tr>", end=' ')
    print('</table>')
#    print '</span></a>'


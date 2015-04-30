# Create your views here.

from django.shortcuts import render
from django.conf import settings
from django.http import HttpResponse, HttpResponseRedirect
from async.forms import ModelForm
import tempfile, os, json, re
# For charts
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.patches import Circle
from itertools import chain
import numpy as np
import json
import mm   #Marmir excel spreadsheet creator

from xpert_su.settings import BASE_DIR

from helper import placement


strat_names = ['0. Baseline (Smear, no Xpert)',
               '1. Xpert for smear-positive',
               '2. Xpert for HIV+',
               '3. Xpert for previously treated',
               '4. Xpert for sm-neg HIV+ or prev tx',
               '5. Xpert for all HIV+ or prev tx',
               '6. Xpert for smear-negative',
               '7. Xpert for all',
               '8. Xpert for all, same-day']

#Top is devel, bottom is production
path_to_xpert_su = BASE_DIR

def splash_page(request):
    return render(request, "splash.html", {})

def home_page(request):
    form = ModelForm()
    return render(request,"index.html", {'form':form})

def about_page(request):
    return render(request,"about.html", {})

def help_page(request):
    return render(request,"help.html", {})

#def homebrew_page(request):
#    return render(request,"homebrew.html", {})

def homebrew_model_page(request):
    full_model = False

    if request.method == 'GET':
        try:
            if request.GET['json_object']:
                full_model = True
        except:
            pass
        
        if full_model:
            #Create the temp file and write the json data to it.
            params = {}
            try:
                params = json.loads(request.GET['json_object'])
            except:
                return render(request,"bad.html", { 'type':'malformed .json' })

            #print params

            params['ip'] = request.META['REMOTE_ADDR'] if request.META['REMOTE_ADDR'] else 'Unknown IP'
      
            params['homebrew'] = True
       
            id, tmp_p = tempfile.mkstemp(suffix='.json',prefix='xpert_hb',
                        dir=settings.MEDIA_ROOT)

            with open(tmp_p, 'w') as fp:
                json.dump(params, fp)

            print params

            #Run the model with the temp file as its first argument
            s = "echo python %s/async/homebrew_model.py " % (path_to_xpert_su)
            s += "%s | at now" % (tmp_p)
            os.system(s)

            change_values = [ ['hivnnew','hiv-new'],
                              ['hivpnew','hiv+new'],
                              ['hivnret','hiv-ret'],
                              ['hivpret','hiv+ret'] ]

            jsonobj = json.loads(request.GET['json_object'])
            
            for v in change_values:
                jsonobj['diag'][v[0]] = jsonobj['diag'][v[1]]
                jsonobj['diag'][v[0]][2] *= 100
                jsonobj['diag'][v[0]][3] *= 100
                #del jsonobj.diag[v[1]]
                jsonobj['dst'][v[0]] = jsonobj['dst'][v[1]]
                #del jsonobj.dst[v[1]]

            jsonobj['ud_dst_tests'] = {}
            for v in [x for x in jsonobj['ud_tests'] if x[:4] == 'DST_']:
                jsonobj['ud_dst_tests'][v[4:]] = jsonobj['ud_tests'][v]
                del jsonobj['ud_tests'][v] 

            #Send the page the file location and the number of strats (always all in this case)
            return render(request,"model.html", {'web_path':"{}{}".format(settings.MEDIA_URL,tmp_p.split("/")[-1]), 'int_select': 9, 'jsonobj' : jsonobj})

    return render(request,"bad.html", { 'type':'homebrew model' })

  
def model_page(request):
    if request.method == "GET":

        full_model = True

        try:
            if request.GET['country'] != 'NAC':
                full_model = False
        except:
            pass 

        if full_model:
            #Run Model
            form = ModelForm(request.GET)
            if form.is_valid():
                strategy = ''
                cd = form.cleaned_data
 
                id, tmp_p = tempfile.mkstemp(suffix='.json',prefix='xpert',
                            dir=settings.MEDIA_ROOT)
                if cd['type_select'] == '0': #if Single Stategy is clicked
                    strategy = cd['int_select']; #run the selected strategy
                else: #'All' is clicked; run all strats.
                    strategy = '9';

                params = {}

                params['int_select'] = int(strategy)
                params['target_inc'] = cd['t_inc']
                params['target_mdr'] = cd['t_mdr']
                params['target_hiv'] = cd['t_hiv']
                params['drug1_cost'] = cd['t_drug1_cost']
                params['drug2_cost'] = cd['t_drug2_cost']
                params['drug3_cost'] = cd['t_drug3_cost']
                params['outpt_cost'] = cd['t_outpt_cost']
                params['sm_cost'] = cd['t_sm_cost']
                params['gxp_cost'] = cd['t_gxp_cost']
                params['sdgxp_cost'] = cd['t_sdgxp_cost']

                params['homebrew'] = False
                params['filename'] = tmp_p

                with open(tmp_p,'w') as fp:
                    json.dump(params, fp)


                s = "echo %s/async/xpert_bg_inter.py " % (path_to_xpert_su)
                s += "%s | at now" % (tmp_p)
                os.system(s)
                return render(request, "model.html", {'web_path':"{}{}".format(settings.MEDIA_URL,tmp_p.split("/")[-1]),
                                                      'int_select':params['int_select']})

        else:
           
            return render(request, "preset.html", {})

    return render(request,"bad.html", { 'type':'model' })

def check_filename (filename):
    status = True

    if len(filename) != 16:
        status = False

    elif filename[:5] != 'xpert':
        status = False

    elif filename[-5:] != '.json':
        status = False

    return status


def excel_sheet(request):
    if request.method == 'GET':

        try:
            filename = request.GET['filename']
        except:
            return render(request, "bad.html", {'type':'no filename in excel GET'})
        
        if not check_filename(filename):
            return render(request, "bad.html", {'type':'bad filename'})

        try:
            with open('%s/media/%s' % (path_to_xpert_su, filename), "r") as fp:
                data = json.load(fp)
        except:
            return render(request, "bad.html", {'type':'no such JSON file'})

        data['intermed'].insert(0, 'Intervention/Scenario')
        
        sheets = ['year1','year2','year3','year4','year5']

        config = { 'header_style':'', 'row_styles': () }

        [ data [str(x)][y].insert (0, strat_names[x]) for x in range (9) for y in sheets ] #add intervention names in

        #total sheet
        total_data = []
        for x in range(9):
            total_data.append([])
            total_data[x].insert(0, strat_names[x])
            for y in range(1,len(data['intermed'])):
                total_data[x].append(0.0)
                for z in sheets:
                    total_data[x][y] += data[str(x)][z][y]
        
        total_sheet = mm.Document (total_data, order=data['intermed'],config_dict=config)

        total_sheet.set_name('Totals')

        mmdocs = [] #This next section could be all list comprehensives but I'm leaving them as for statements
                    #to reduce visual complexity.  Also programmer complexity.

        for x in range(len(sheets)):
            mmdocs.append(mm.Document ([data[str(y)][sheets[x]] for y in range(9) ],order=data['intermed'],config_dict=config))
            mmdocs[x].set_name(sheets[x])

        #Add all excel sheets to the first sheet
        [ total_sheet.add_child(mmdocs[x]) for x in range(len(sheets)) ]

        #write to a file:
        write_filename = '/media/excel_%s.xls' % ( filename.split('.json')[0] )
        total_sheet.write('%s%s' % (path_to_xpert_su, write_filename) )

        #Let Apache serve it
        return HttpResponseRedirect ( write_filename )

    else:
        return render(request, "bad.html", {'type':'excel'})

def dgraph1(request):
    p_val = request.GET

    try:
        cmore = [float(p_val['c{}'.format(x)]) for x in range(9)]
        inc = [float(p_val['i{}'.format(x)]) for x in range(9)]

    except:
        return render(request,"bad.html", { 'type':'dot graph 1' })
       
    which = -1
    try:
        which = int(p_val['ref'])
    except:
        pass

    s8name = strat_names[8]
    try:
        s8name = p_val['s8name']
    except:
        pass
                
    fig = Figure(facecolor='white', dpi=150)
    canvas = FigureCanvas(fig)

    ax1 = fig.add_subplot(111)
    
    ax1.set_xlabel('% decrease in TB incidence at year 5 (more effective -->)',{'fontweight':'bold'})
    ax1.set_ylabel('% increase in cost at year 5 (more costly -->)',{'fontweight':'bold'})

    points = [ ax1.plot(inc[i],cmore[i],color='blue',marker='.',linestyle='None') for i in xrange(9) ]

    ax1.set_title("Change (%) in cost and TB incidence at year 5",{'fontweight':'bold'})

    ax1.axhline(color='r')
    ax1.axvline(color='r')

#    ax1.plot ([0],[0],marker='*',c='black',linestyle='None')    

    #Add annotations, simple collision detection
    point_list = []
    for i in xrange(9):
        dx,dy = placement(inc[i],cmore[i],point_list,ax1)
        point_list.append([inc[i],cmore[i],dx,dy])
        ref_color = '#000000'
        if i == which:
            ref_color = 'red'
        ax1.annotate(str(i), xy=(inc[i],cmore[i]), 
                     xytext=(inc[i]+dx,cmore[i]+dy),
                     color=ref_color,
                     arrowprops=dict(color='red',arrowstyle="->"))

    fig.tight_layout()
    
    ax1.legend ((points[0][0],points[1][0],points[2][0],points[3][0],points[4][0],points[5][0],points[6][0],points[7][0],points[8][0]), 
                (strat_names[0],
                 strat_names[1],
                 strat_names[2],
                 strat_names[3],
                 strat_names[4],
                 strat_names[5],
                 strat_names[6],
                 strat_names[7],
                 "8. {}".format(s8name)),
                 loc='best', numpoints=1,prop={'size':'9'})

    response=HttpResponse(content_type='image/png')
    canvas.print_png(response,facecolor=fig.get_facecolor())
    return response

def dgraph2(request):
    p_val = request.GET

    try:
        cmore = [float(p_val['c{}'.format(x)]) for x in range(9)]
        mdr = [float(p_val['m{}'.format(x)]) for x in range(9)]
    except:
        return render(request,"bad.html",{ 'type':'dot graph 2' })
           
    which = -1
    try:
        which = int(p_val['ref'])
    except:
        pass

    s8name = strat_names[8]
    try:
        s8name = p_val['s8name']
    except:
        pass

    fig = Figure(facecolor='white', dpi=150)
    canvas = FigureCanvas(fig)

    ax2 = fig.add_subplot(111)

    ax2.set_xlabel('% decrease in MDR-TB incidence at year 5 (more effective -->)',{'fontweight':'bold'})
    ax2.set_ylabel('% increase in cost at year 5 (more costly -->)',{'fontweight':'bold'})

    points = [ ax2.plot(mdr[i],cmore[i],color='blue',marker='.',linestyle='None') for i in xrange(9) ]

    ax2.set_title("Change (%) in cost and MDR incidence at year 5",{'fontweight':'bold'})

    ax2.axhline(color='r')
    ax2.axvline(color='r')

    ax2.plot ([0],[0],marker='*',c='black',linestyle='None')


    point_list = []
    for i in xrange(9):
        dx,dy = placement(mdr[i],cmore[i],point_list,ax2)
        point_list.append([mdr[i],cmore[i],dx,dy])
        ref_color = '#000000'
        if i == which:
            ref_color = 'red'
        ax2.annotate(str(i), xy=(mdr[i],cmore[i]),
                     xytext=(mdr[i]+dx,cmore[i]+dy),
                     color=ref_color,
                     arrowprops=dict(color='red',arrowstyle="->"))    

    fig.tight_layout()

    ax2.legend ((points[0][0],points[1][0],points[2][0],points[3][0],points[4][0],points[5][0],points[6][0],points[7][0],points[8][0]), 
                (strat_names[0],
                 strat_names[1],
                 strat_names[2],
                 strat_names[3],
                 strat_names[4],
                 strat_names[5],
                 strat_names[6],
                 strat_names[7],
                 "8. {}".format(s8name)),
                 loc='best', numpoints=1,prop={'size':'9'})

    response=HttpResponse(content_type='image/png')
    canvas.print_png(response,facecolor=fig.get_facecolor())
    return response

def dgraph3(request):
    p_val = request.GET
    
    try:
        cmore = [float(p_val['cb{}'.format(x)]) for x in xrange(1,9)]
        inc   = [float(p_val['tb{}'.format(x)]) for x in xrange(1,9)]
        c_hi  = [float(p_val['cbh{}'.format(x)]) for x in xrange(1,9)]
        c_lo  = [float(p_val['cbl{}'.format(x)]) for x in xrange(1,9)]
        inc_hi= [float(p_val['tbh{}'.format(x)]) for x in xrange(1,9)]
        inc_lo= [float(p_val['tbl{}'.format(x)]) for x in xrange(1,9)]
    except:
        return render(request,"bad.html", { 'type':'dot graph 1' })
   
    which = -1
    try:
        which = int(p_val['ref'])
    except:
        pass
    
    if p_val['flag'] == 'E':
        try: 
            e2i = [float(p_val['e1%dt'%(x)]) for x in xrange(1,9)]
            e2c = [float(p_val['e1%dc'%(x)]) for x in xrange(1,9)]
        except:
            return render(request,"bad.html", { 'type':'uncert_e bargraph' })
    elif p_val['flag'] == 'P':
        try:
            p2i = [float(p_val['p1%dt'%(x)]) for x in xrange(1,9)]
            p2c = [float(p_val['p1%dc'%(x)]) for x in xrange(1,9)]
        except:
            return render(request,"bad.html", { 'type':'uncert_p bargraph' })
    elif p_val['flag'] == 'R':
        try:
            r2i = [float(p_val['r1%dt'%(x)]) for x in xrange(1,9)]
            r2c = [float(p_val['r1%dc'%(x)]) for x in xrange(1,9)]
        except:
            return render(request,"bad.html", { 'type':'uncert_r bargraph' })

 
    fig = Figure(facecolor='white', dpi=150) #
    canvas = FigureCanvas(fig)

    ax = fig.add_subplot(111)

    ax.set_xlabel('% decrease in TB incidence at year 5 (more effective -->)',{'fontweight':'bold'})
    ax.set_ylabel('% increase in cost at year 5 (more costly -->)',{'fontweight':'bold'})
    ax.axvline(color='r')
    
    #n_colours = ['#fb9902','#fd5308','#fe2712','#a7194b','#8601af','#3d01a4','#0247fe','#0392ce','#66b032']

    n_colours = ['#000000','#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#a65628','#f781bf','#999999']

    colours = ['#CC0000','#EC6002','#CC9999','#FFFF99','#CCFF00','#02C5F4','#0188A9','#006699','#2E014B']

    t_delta_x = 0.25
    t_delta_y = 2

    ax.text(0,0,'0')

    for i in xrange(8):
        ax.text(inc[i]+t_delta_x,cmore[i]+t_delta_y,str(i+1))

    delta = 0.00002
    
    errorbars = []
    
    errorbars.append(ax.plot([0],[0],color='black',marker='*'))

    for i in xrange(8):
        centre_color = n_colours[i+1]
        if which-2 == i:
            centre_color = 'black'
        errorbars.append(ax.errorbar([inc[i]],[cmore[i]],
                  xerr=[[inc[i] - inc_lo[i]],[inc_hi[i] - inc[i]]],
                  yerr=[[cmore[i] - c_lo[i]],[c_hi[i] - cmore[i]]],
                  color=centre_color, #n_colours[i+1]
                  ecolor=centre_color,
                  linewidth=2))
    ax.axhline(color='r')

    ax.plot ([0],[0],marker='*',c='black',linestyle='None')

    if p_val['flag'] == 'E':
        for x in range(8):
            centre_color = n_colours[x+1]
            if which-2 == x:
                centre_color = 'black'
            ax.plot([e2i[x]],[e2c[x]],marker='o',c=centre_color,linestyle='None')
    elif p_val['flag'] == 'P':
        for x in range(8):
            centre_color = n_colours[x+1]
            if which-2 == x:
                centre_color = 'black'
            ax.plot([p2i[x]],[p2c[x]],marker='o',c=centre_color,linestyle='None')
    elif p_val['flag'] == 'R':
        for x in range(8):
            centre_color = n_colours[x+1]
            if which-2 == x:
                centre_color = 'black'
            ax.plot([r2i[x]],[r2c[x]],marker='o',c=centre_color,linestyle='None')

    ax.set_title("Change (%) in cost and TB incidence at year 5",{'fontweight':'bold'})

    box = ax.get_position()
    ax.set_position([box.x0, box.y0-0.05, box.width, box.height+0.05])

    #ax.legend(errorbars,range(1,9),numpoints=1,bbox_to_anchor=(0.98,-0.08), ncol=9,columnspacing=1)

    fig.tight_layout()
  
    ax.legend ((errorbars[0][0],errorbars[1][0],errorbars[2][0],errorbars[3][0],errorbars[4][0],errorbars[5][0],errorbars[6][0],errorbars[7][0],errorbars[8][0]), 
               (strat_names[0],
                strat_names[1],
                strat_names[2],
                strat_names[3],
                strat_names[4],
                strat_names[5],
                strat_names[6],
                strat_names[7],
                strat_names[8]),
                loc='best', numpoints=1,prop={'size':'9'})

    response=HttpResponse(content_type='image/png')
    canvas.print_png(response,facecolor=fig.get_facecolor())
    return response

def dgraph4(request):
    p_val = request.GET

    try:
        cmore = [float(p_val['cb{}'.format(x)]) for x in range(1,9)]
        mdr = [float(p_val['mb{}'.format(x)]) for x in range(1,9)]
        c_hi  = [float(p_val['cbh{}'.format(x)]) for x in xrange(1,9)]
        c_lo  = [float(p_val['cbl{}'.format(x)]) for x in xrange(1,9)]
        mdr_hi= [float(p_val['mbh{}'.format(x)]) for x in xrange(1,9)]
        mdr_lo= [float(p_val['mbl{}'.format(x)]) for x in xrange(1,9)]

    except:
        return render(request,"bad.html",{ 'type':'dot graph 2' })

    which = -1
    try:
        which = int(p_val['ref'])
    except:
        pass
   
    if p_val['flag'] == 'E':
        try: 
            e2m = [float(p_val['e1%dm'%(x)]) for x in xrange(1,9)]
            e2c = [float(p_val['e1%dc'%(x)]) for x in xrange(1,9)]
        except:
            return render(request,"bad.html", { 'type':'uncert_e bargraph' })
    elif p_val['flag'] == 'P':
        try:
            p2m = [float(p_val['p1%dm'%(x)]) for x in xrange(1,9)]
            p2c = [float(p_val['p1%dc'%(x)]) for x in xrange(1,9)]
        except:
            return render(request,"bad.html", { 'type':'uncert_p bargraph' })
    elif p_val['flag'] == 'R':
        try:
            r2m = [float(p_val['r1%dm'%(x)]) for x in xrange(1,9)]
            r2c = [float(p_val['r1%dc'%(x)]) for x in xrange(1,9)]
        except:
            return render(request,"bad.html", { 'type':'uncert_r bargraph' })

    fig = Figure(facecolor='white', dpi=150)
    canvas = FigureCanvas(fig)

    ax = fig.add_subplot(111)

    ax.set_xlabel('% decrease in MDR-TB incidence at year 5 (more effective -->)',{'fontweight':'bold'})
    ax.set_ylabel('% increase in cost at year 5 (more costly -->)',{'fontweight':'bold'})
    ax.axvline(color='r')

    #n_colours = ['#fb9902','#fd5308','#fe2712','#a7194b','#8601af','#3d01a4','#0247fe','#0392ce','#66b032']

    n_colours = ['#000000','#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#a65628','#f781bf','#999999']

    colours = ['#CC0000','#EC6002','#CC9999','#FFFF99','#CCFF00','#02C5F4','#0188A9','#006699','#2E014B']

    t_delta_y = 2
    t_delta_x = 0.25

    ax.text(0,0,'0')

    for i in xrange(8):
        ax.text(mdr[i]+t_delta_x,cmore[i]+t_delta_y,str(i+1))

    delta = 0.00002

    errorbars = []

    errorbars.append(ax.plot([0],[0],color='black', marker='*'));
    
    for i in xrange(8):
        centre_color = n_colours[i+1]
        if which-2 == i:
            centre_color = 'black'
        errorbars.append(ax.errorbar([mdr[i]],[cmore[i]],
                  xerr=[[mdr[i] - mdr_lo[i]],[mdr_hi[i] - mdr[i]]],
                  yerr=[[cmore[i] - c_lo[i]],[c_hi[i] - cmore[i]]],
                  color=centre_color,
                  ecolor=centre_color,
                  linewidth=2))
    ax.axhline(color='r') 
   
    ax.plot ([0],[0],marker='*',c='black',linestyle='None')    

    if p_val['flag'] == 'E':
        for x in range(8):
            ax.plot([e2m[x]],[e2c[x]],marker='o',c=n_colours[x+1],linestyle='None')
    elif p_val['flag'] == 'P':
        for x in range(8):
            ax.plot([p2m[x]],[p2c[x]],marker='o',c=n_colours[x+1],linestyle='None')
    elif p_val['flag'] == 'R':
        for x in range(8):
            ax.plot([r2m[x]],[r2c[x]],marker='o',c=n_colours[x+1],linestyle='None')

    ax.set_title("Change (%) in cost and MDR incidence at year 5",{'fontweight':'bold'})
 
    box = ax.get_position()
    ax.set_position([box.x0, box.y0-0.05, box.width, box.height+0.05])

    #ax.legend(errorbars,range(1,9),numpoints=1,bbox_to_anchor=(0.98,-0.08),ncol=9,columnspacing=1)
    fig.tight_layout()

    ax.legend ((errorbars[0][0],errorbars[1][0],errorbars[2][0],errorbars[3][0],errorbars[4][0],errorbars[5][0],errorbars[6][0],errorbars[7][0],errorbars[8][0]), 
               (strat_names[0],
                strat_names[1],
                strat_names[2],
                strat_names[3],
                strat_names[4],
                strat_names[5],
                strat_names[6],
                strat_names[7],
                strat_names[8]),
                loc='best', numpoints=1,prop={'size':'9'})

    response=HttpResponse(content_type='image/png')
    canvas.print_png(response,facecolor=fig.get_facecolor())
    return response

def bargraph(request):
    p = request.GET    

    try:
        d = [(float(p['d10']), float(p['d11']), float(p['d12']), float(p['d13']), float(p['d14'])),
             (float(p['d20']), float(p['d21']), float(p['d22']), float(p['d23']), float(p['d24'])),
             (float(p['d30']), float(p['d31']), float(p['d32']), float(p['d33']), float(p['d34'])),
             (float(p['d40']), float(p['d41']), float(p['d42']), float(p['d43']), float(p['d44'])),
             (float(p['d50']), float(p['d51']), float(p['d52']), float(p['d53']), float(p['d54'])),
             (float(p['d60']), float(p['d61']), float(p['d62']), float(p['d63']), float(p['d64'])),
             (float(p['d70']), float(p['d71']), float(p['d72']), float(p['d73']), float(p['d74'])),
             (float(p['d80']), float(p['d81']), float(p['d82']), float(p['d83']), float(p['d84']))]
    except:
        return render(request,"bad.html", { 'type':'bargraph' })
    
    s8name = strat_names[8]
    try:
        s8name = p['s8name']
    except:
        pass

    colors = ["grey","blue","green","yellow","red"]

    cdata = []

    for array in zip(*d): #Invert the numbers for proper display
        cdata.append([-1 * x for x in array])

    loc = np.arange(len(cdata[0]))
 
    width = 0.15

    fig = Figure(facecolor='white', dpi=150, figsize=(8,6))
    canvas = FigureCanvas(fig)

    ax = fig.add_subplot(111)

    rect = [ax.bar(loc+width*i, cdata[i], width, color=colors[i]) 
            for i in range(len(cdata))]
    
    max_height = max(chain(cdata[0],cdata[1],cdata[2],cdata[3],cdata[4])) + 10
    min_height = min(chain(cdata[0],cdata[1],cdata[2],cdata[3],cdata[4])) - 10

    print max_height, min_height

    ax.set_ylim(min_height,max_height)
    ax.set_xlim(-width*4, len(loc) +(4*width))

    ax.set_xticks(loc + (2.5*width))

    mod_strat_names = strat_names[1:-1] + ["8. {}".format(s8name)]

    ax.set_xticklabels(mod_strat_names, rotation='30', size='small', stretch='condensed',
                       ha='right' )

    ax.legend ((rect[0][0], rect[1][0], rect[2][0], rect[3][0], rect[4][0]),
                ("TB Incidence", "MDR Incidence", "TB Mortality","Cost at Year 1", "Cost at Year 5"),loc='best',
                prop={'size':'8'})

    ax.set_title ("Side-by-Side Comparison of the Strategies Using 5 Key Metrics", {'fontweight':'bold'})
    ax.axhline(color='black')

    ax.set_ylabel('percentage change from baseline',{'fontweight':'bold'})

    fig.tight_layout()

    response=HttpResponse(content_type='image/png')
    canvas.print_png(response,facecolor=fig.get_facecolor())
    return response

def bargraph_uncert(request):
    p = request.GET    

    #BASIC TEST

    try:
        ca  = [float(p['ca{}'.format(x)]) for x in range(1,9)]
        cal = [float(p['cal{}'.format(x)]) for x in range(1,9)]
        cah = [float(p['cah{}'.format(x)]) for x in range(1,9)]
        cb  = [float(p['cb{}'.format(x)]) for x in range(1,9)]
        cbl = [float(p['cbl{}'.format(x)]) for x in range(1,9)]
        cbh = [float(p['cbh{}'.format(x)]) for x in range(1,9)]
        tb  = [float(p['tb{}'.format(x)]) for x in xrange(1,9)]
        tbl = [float(p['tbl{}'.format(x)]) for x in xrange(1,9)]
        tbh = [float(p['tbh{}'.format(x)]) for x in xrange(1,9)]
        mb  = [float(p['mb{}'.format(x)]) for x in xrange(1,9)]
        mbl = [float(p['mbl{}'.format(x)]) for x in xrange(1,9)]
        mbh = [float(p['mbh{}'.format(x)]) for x in range(1,9)]
        ob  = [float(p['ob{}'.format(x)]) for x in range(1,9)]
        obl = [float(p['obl{}'.format(x)]) for x in range(1,9)]
        obh = [float(p['obh{}'.format(x)]) for x in range(1,9)]
    except:
        return render(request,"bad.html", { 'type':'bargraph' })

    colors = ["grey","blue","green","yellow","red"]
    
    pdata = [ tb, mb, ob, ca, cb ]
    
    #ndata = zip(*d)

    height = max(chain(tbh,mbh,obh,cah,cbh,tb,mb,ob,ca,cb)) + 10

    #loc = np.arange(len(ndata[0])) 
    loc = np.arange(len(ca))

    width = 0.15

    fig = Figure(facecolor='white', dpi=150, figsize=(8,6))
    canvas = FigureCanvas(fig)

    ax = fig.add_subplot(111)

    rect = [ax.bar(loc+width*i, pdata[i], width, color=colors[i]) 
            for i in range(len(pdata))]

    [ax.errorbar([i+(width*.5)], tb[i], yerr=[[tb[i] - tbl[i]],[tbh[i] - tb[i]]],
                   color='black',
                   ecolor='black',
                   linewidth=1)
                   for i in xrange(8)] #TBIncidence Errorbar
    [ax.errorbar([i+(width*1.5)], mb[i], yerr=[[mb[i] - mbl[i]],[mbh[i] - mb[i]]],
                   color='black',
                   ecolor='black',
                   linewidth=1)
                   for i in xrange(8)] #MDRIncidence Errorbar
    [ax.errorbar([i+(width*2.5)], ob[i], yerr=[[ob[i] - obl[i]],[obh[i] - ob[i]]],
                   color='black',
                   ecolor='black',
                   linewidth=1)
                   for i in xrange(8)] #Mortality Errorbar
    [ax.errorbar([i+(width*3.5)], ca[i], yerr=[[ca[i] - cal[i]],[cah[i] - ca[i]]],
                   color='black',
                   ecolor='black',
                   linewidth=1)
                   for i in xrange(8)] #Cost Yr 1 Errorbar
    [ax.errorbar([i+(width*4.5)], cb[i], yerr=[[cb[i] - cbl[i]],[cbh[i] - cb[i]]],
                   color='black',
                   ecolor='black',
                   linewidth=1)
                   for i in xrange(8)] #Cost Yr 5 Errorbar

    ax.set_ylim(-50,height)
    ax.set_xlim(-width*4, len(loc) +(4*width))

    ax.set_xticks(loc + (2.5*width))

    ax.set_xticklabels(strat_names[1:], rotation='30', size='small', stretch='condensed',
                       ha='right' )

    ax.legend ((rect[0][0], rect[1][0], rect[2][0], rect[3][0], rect[4][0] ),
                ("TB Incidence", "MDR Incidence", "TB Mortality","Cost at Year 1", "Cost at Year 5"),loc='best',
                prop={'size':'8'})

    ax.set_title ("Side-by-Side Comparison of the Strategies Using 5 Key Metrics", {'fontweight':'bold'})
    ax.axhline(color='black')

    ax.set_ylabel('percentage change from baseline',{'fontweight':'bold'})

    fig.tight_layout()

    response=HttpResponse(content_type='image/png')
    canvas.print_png(response,facecolor=fig.get_facecolor())
    return response

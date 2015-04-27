$(document).ready( function () {
    
    var XPT_u = {};

    XPT_u.data = {};
    /*XPT_u.load_mo_value = '';*/

    XPT_u.which = 0;


    var QueryString = function(name) {
        name = name.replace(/[\[]/, "\\[").replace(/[\]]/, "\\]");
        var regex = new RegExp("[\\?&]" + name + "=([^&#]*)"),
            results = regex.exec(location.search);
        return results == null ? "" : decodeURIComponent(results[1].replace(/\+/g, " "));
    };
     
    XPT_u.iso3 = QueryString('country');

    (function (jsonfile) {
        $.ajax ({
            url:jsonfile,
            async:false,
            dataType:'json',
            success: function (data) {
                XPT_u.data = data;
            },
            error: function () {
                console.log('JSON not loaded');
            }

        });
    })('/static/data_Sept3/' + XPT_u.iso3 + '.json');

    (function (jsonfile) {
        $.ajax ({
            url:jsonfile,
            async:false,
            dataType:'json',
            success: function (data) {
                XPT_u.mOver_data = data;
            },
            error: function () {
                console.log('Could not load mouse over data');
            }
        });
    })('/static/mouseOver.json');

    function round_to( val, digits ) {
        return Math.round(val * Math.pow(10,digits)) / Math.pow(10,digits);
    }

    XPT_u.add_brackets = function ( val ) {
        var addline;
        addline = val.toString();
        if (val < 0)
            addline = "(" + addline + ")";
        return addline;
    };

    /* Functions to load data onto page */

    XPT_u.load_mo_value = function ( which ) {

        var point = parseFloat(XPT_u.data[which + '_point']);
        var hi    = parseFloat(XPT_u.data[which + '_hi']);
        var lo    = parseFloat(XPT_u.data[which + '_lo']);
        var delta = 0.00002;

        if ( Math.abs(point - 0.0) < delta || isNaN(point)) {
            $("[name='" + which + "_point']").text('Missing');
        } else {
            var txt = '';
            if (which === 'tb') {
                txt = ' per 100,000';
            }
            $("[name='" + which + "_point']").text(round_to(point, 1) + txt);
        }

        if ( Math.abs(lo - 0.0) < delta || isNaN(lo)) {
            $("[name='" + which + "_lo']").text('Missing');
        } else {
            $("[name='" + which + "_lo']").text(round_to(lo, 1));
        }
      
        if ( Math.abs(hi - 0.0) < delta || isNaN(hi)) {
            $("[name='" + which + "_hi']").text('Missing');
        } else {
            $("[name='" + which + "_hi']").text(round_to(hi, 1));
        }

    };

    XPT_u.load_co_value = function () {
        var dict = {};
        dict.drg_cst = parseFloat(XPT_u.data['drug1_cost']); /* 500 */
        dict.drg2_cst = dict.drg_cst * 2;
        dict.drg3_cst = dict.drg_cst * 10;
        dict.opt_cst = parseFloat(XPT_u.data['outpt_cost']); /* 10 */
        dict.sm_cst = (dict.opt_cst / 5.0) > 2.0?(dict.opt_cst / 5.0):2.0;
        dict.sd_cst = dict.sm_cst * 5;
        dict.gxp_cst = 15;
        dict.sdgxp_cst = (dict.sd_cst - dict.sm_cst) + dict.gxp_cst;
        dict.cx_cst = dict.sm_cst * 10;
        dict.dst_cst = dict.sm_cst * 20;
        dict.mods_cst = dict.sm_cst * 2.5;
        return dict;
    };

    XPT_u.load_name = function () {
        $("[name='preset_fullname']").text(XPT_u.data['fullname']);
        $("[name='preset_gdp']").text(round_to(XPT_u.data['gdp'],2));
    };

    XPT_u.load_m_options = function () {
        /* TB incidence */
        XPT_u.load_mo_value('tb');
        XPT_u.load_mo_value('aids');
        XPT_u.load_mo_value('mdr');
    };

    XPT_u.load_c_options = function () {
        var values = XPT_u.load_co_value();
        var name_array = ['drg_cst','drg2_cst','drg3_cst','opt_cst',
                          'sm_cst','gxp_cst','sdgxp_cst'];
        $.each(name_array, function (v, d) {
            $("[name='" + d + "']").text("$"+round_to(parseFloat(values[d]),2).toString());
        });
    };

    XPT_u.load_graph_var = function (name, xdigit, neg) {
        var store = [[],[],[]];
        var letter = (name[name.length - 1] === '1')?'a':'b';
        for (var i = 1; i < 9; i ++ ) {
            if (neg) {
                store[0].push(xdigit+letter+i.toString()+"="+(parseFloat(XPT_u.data[i.toString()][name]['value'])*-1).toString()+'&');
                store[1].push(xdigit+letter+'l'+i.toString()+"="+(parseFloat(XPT_u.data[i.toString()][name]['lo'])*-1).toString()+'&');
                store[2].push(xdigit+letter+'h'+i.toString()+"="+(parseFloat(XPT_u.data[i.toString()][name]['hi'])*-1).toString()+'&');
            } else {
                store[0].push(xdigit+letter+i.toString()+"="+XPT_u.data[i.toString()][name]['value']+'&');
                store[1].push(xdigit+letter+'l'+i.toString()+"="+XPT_u.data[i.toString()][name]['lo']+'&');
                store[2].push(xdigit+letter+'h'+i.toString()+"="+XPT_u.data[i.toString()][name]['hi']+'&');
            }
        }
        return store;

    };

    XPT_u.load_core_var = function (name, xdigit, inc) {
        var store = [];
        for (var i = 1; i < 9; i ++) { 
            store.push(name+i.toString()+inc[0]+"="+XPT_u.data[i.toString()][name][inc]+'&');
            store.push(name+i.toString()+'c='+XPT_u.data[i.toString()][name]['cost']+'&');
        }
        return store;
    };

    XPT_u.create_graph = function (path, xvar, yvar, div_id, flag) {
        var fullURL = '';

        var xarr = XPT_u.load_graph_var(xvar, xvar[0]);

        var yarr = XPT_u.load_graph_var(yvar, yvar[0]);
        
        var supp = '';
        
        var descrip = '#' + div_id.split('_')[0] + '_type';
        var radio_selected;

        switch (flag) {
            case 'B':
                $(descrip).text('');
                break;
            case 'E':
                var e_core = XPT_u.load_core_var ('e1', 'e', yvar);
                supp = e_core.join('');
                $(descrip).text('Dots are doubling of empiric treatment');
                break;
            case 'P':
                var p_core = XPT_u.load_core_var ('p1', 'p', yvar);
                supp = p_core.join('');
                $(descrip).text('Dots are doubling of pre-diagnostic delay');
                break;
            case 'R':
                var r_core = XPT_u.load_core_var ('r1', 'r', yvar);
                supp = r_core.join('');
                $(descrip).text('Dots are doubling of reactivation');
                break;
            default:
                break;
        }
        
        radio_selected = $('[name="rbg_' + yvar + '"]:checked').val();

        fullURL = path + xarr[0].join('') +
                         xarr[1].join('') +
                         xarr[2].join('') +
                         yarr[0].join('') +
                         yarr[1].join('') +
                         yarr[2].join('') +
                         supp + "flag=" + flag +
                         "&ref=" + (parseInt(radio_selected,10)+1).toString();

        $('#'+div_id+' img').attr('src',fullURL);
        
        $('#graph'+div_id.split('sp')[1]+"_frame").find('a').each(function (index, data) { $(data).remove(); });

        $('#sp'+div_id.split('sp')[1]+"_table").before('<a href="' + fullURL + '" class="btn btn-sm btn-primary d_plot_button" download="' + xvar + yvar + '.png">Download Plot</a>');
    };

    XPT_u.create_table = function ( type, div_id, flag, order_by ) {

        var names = ['Xpert for smear+',
                   'Xpert for HIV+',
                   'Xpert for prev-tx',
                   'Xpert for sm-neg HIV+ or prev-tx',
                   'Xpert for all HIV+ or prev-tx',
                   'Xpert for smear negative',
                   'Xpert for all',
                   'Xpert for all, same-day'
                  ];

/*        var n_colours = ['#fb9902','#fd5308','#fe2712','#a7194b','#8601af',
                         '#3d01a4','#0247fe','#0392ce','#66b032'];*/
        var n_colours = ['#000000','#e41a1c','#377eb8','#4daf4a','#984ea3',
                         '#ff7f00','#a65628','#f781bf','#999999'];

        var colours = ['#CC0000','#EC6002','#CC9999','#FFFF99','#CCFF00',
                       '#02C5F4','#0188A9','#006699','#2E014B'];
        var t_str = '';
        var addline;
    
        var moveHeight = [ 0,55,0,0,55,55,55,75,55 ];

        var row_colour;

        var baseline = {'cost_val':0,'cost_lo':0,'cost_hi':0,'type_val':0,'type_lo':0,'type_hi':0};

        if (order_by != 0) { 
            baseline['cost_val'] = XPT_u.data[order_by.toString()]['cost5']['value'];
            baseline['cost_lo']  = XPT_u.data[order_by.toString()]['cost5']['lo'];
            baseline['cost_hi']  = XPT_u.data[order_by.toString()]['cost5']['hi'];
            baseline['type_val'] = XPT_u.data[order_by.toString()][type]['value'];
            baseline['type_lo']  = XPT_u.data[order_by.toString()][type]['lo'];
            baseline['type_hi']  = XPT_u.data[order_by.toString()][type]['hi'];
        }

        t_str += "<table class='leg_disp'>";

        if (type === 'tbi') {
            t_str += "<tr><th colspan='6'>Percent change in cost and TB incidence compared to Baseline (Smear) at year 5</th>";
            t_str += "<th>Ref. Std.</td></tr>";
            t_str += "<tr><th>Number</th><th>Name</th><th>Chg in Cost</th><th>95% Uncert.</th><th>Chg in TB Incidence</th><th>95% Uncert.</th>";
            t_str += "<th><input id='rb_tbi_0' type='radio' name='rgb_tbi' value='0'";
        } else {
            t_str += "<tr><th colspan='6'>Percent change in cost and MDR incidence compared to Baseline (Smear) at year 5</th>";
            t_str += "<th>Ref. Std.</th></tr>";
            t_str += "<tr><th>Number</th><th>Name</th><th>Chg in Cost</th><th>95% Uncert.</th><th>Chg in MDR Incidence</th><th>95% Uncert.</th>";
            t_str += "<th><input id='rb_mdr_0' type='radio' name='rgb_mdr' value='0'";
        }

        if (order_by == 0) {
            t_str += " checked='checked'";
        }
        t_str += "></th></tr>";

        for (var x = 1; x < 9; x ++) {
            row_colour = n_colours[x]

            if (order_by == x) {
                row_colour = '#000000';
            }
            t_str += "<tr id='td_" + x + "' class='tooltipsrc' style='color:"+row_colour+";'><th class='text-center' style='color:"+row_colour+";'>"+(x).toString()+"</th><td>"+names[x-1]+"</td>";

            t_str += "<td>"+round_to(XPT_u.data[x.toString()]['cost5']['value']-baseline['cost_val'],0).toString()+"%</td>";
            
            addline = XPT_u.add_brackets(round_to(XPT_u.data[x.toString()]['cost5']['lo']-baseline['cost_lo'],0));
            t_str += "<td>"+ addline +"% - ";

            addline = XPT_u.add_brackets(round_to(XPT_u.data[x.toString()]['cost5']['hi']-baseline['cost_hi'],0));
            t_str += addline + "%</td>";

            t_str += "<td>"+round_to(XPT_u.data[x.toString()][type]['value']-baseline['type_val'],0).toString()+"%</td>";

            
            addline = XPT_u.add_brackets(round_to(XPT_u.data[x.toString()][type]['lo']-baseline['type_lo'],0));
            t_str += "<td>"+ addline +"% - ";

            addline = XPT_u.add_brackets(round_to(XPT_u.data[x.toString()][type]['hi']-baseline['type_hi'],0));

            t_str += addline +"%</td>";
            t_str += "<td class='radio_buttons'><input id='rb_" + type + "_" + x + "' type='radio' name='rbg_"+type+"' value='" + x + "'";
            if (order_by == x) {
                t_str += " checked='checked'";
            }
            t_str+= "></td></tr>";
        }
        
        t_str += "</table>";
        $('#'+div_id).empty();
        $('#'+div_id).append(t_str);

        XPT_u.load_events();

        if (type =='tbi') {
            XPT_u.load_graph3( flag );
        } else {
            XPT_u.load_graph4( flag );
        }

    };

    XPT_u.create_barbox = function ( path, div_id, flag ) {
        var fullURL = '';

        var mdr = XPT_u.load_graph_var('mdr', 'm', true);
        var inc = XPT_u.load_graph_var('tbi', 't', true);
        var cs5 = XPT_u.load_graph_var('cost5', 'c');
        var cs1 = XPT_u.load_graph_var('cost1', 'c');
        var mor = XPT_u.load_graph_var('mortality','o', true);
    
        fullURL = path + mdr[0].join('') + mdr[1].join('') + mdr[2].join('') +
                         inc[0].join('') + inc[1].join('') + inc[2].join('') +
                         mor[0].join('') + mor[1].join('') + mor[2].join('') +
                         cs5[0].join('') + cs5[1].join('') + cs5[2].join('') +
                         cs1[0].join('') + cs1[1].join('') + cs1[2].join('');

        $('#'+div_id).empty();
        $('#'+div_id).append("<img src='"+fullURL+"' class='pre_graph' alt='Barchart' />");
        $("#"+div_id).prepend("<a href='" +fullURL + "' class='btn btn-sm btn-primary d_plot_button' download='barchart.png'>Download Plot</a>");

    };

    XPT_u.create_comp_table = function ( tb ) {
        var max_g = 0.0;
        var max_r = 0.0;
        var h_a = ['Chg in total inc.',
                   'Chg in MDR inc.',
                   'Chg in TB mort.',
                   'Cost Chg Yr1',
                   'Cost Chg Yr5'];
        var s_a = ['Xpert for smear+',
                   'Xpert for HIV+',
                   'Xpert for prev-tx',
                   'Xpert for sm-neg HIV+ or prev-tx',
                   'Xpert for all HIV+ or prev-tx',
                   'Xpert for smear negative',
                   'Xpert for all',
                   'Xpert for all, same-day.'
                  ];
        var v_names = ['tbi','mdr','mortality','cost1','cost5'];
        var t_string = ""

        t_string+="<table class='table disp'>";
        t_string+="<tr><th>Scenario #</th><th>Name";
        for (var j=0;j<5;j++) {
            t_string +="</th><th>"+h_a[j];
        }
        t_string += "</th></tr>"

        for (var j=0;j<8;j++) {
           
            t_string+="<th>"+(j+2).toString()+"</th><th>"+s_a[j]+"</th>";
            for (var x=0;x<5;x++) {
                var data_value = parseFloat(XPT_u.data[(j+1).toString()][v_names[x]]['value']);
                if (v_names[x] !== 'cost1' && v_names[x] !== 'cost5') {
                    data_value = -1 * data_value;
                }
                if (data_value < 0) {
                    if (Math.abs(data_value) > max_g) {
                        max_g = Math.abs(data_value);
                    }
                    t_string+="<td class='green'>" + round_to(data_value*-1,0).toString() + "%</td>";  
                } else {
                    if (data_value > max_r) {
                        max_r = data_value;
                    }
                    t_string+="<td class='red'>"+ round_to(data_value,0).toString() +"%</td>";
                }
            }
            t_string+="</tr>";
        }                    
        t_string+="<tr><td colspan='7' class='comment'>";
        t_string+="All values indicate pecentage change from baseline scenario (smear)</td></tr>"
        t_string+="<tr><td colspan='7' class='comment'>";
        t_string+="Green values represent reductions, red values represent increases</td></tr>"
         
        /*tb.find("table").append(t_string);*/
        t_string += "</table>";

        tb.append(t_string);

        tb.find("td").each(function () {
            var start = 50; /*96*/
            var max = 255;  /*230*/ /*Changed to provide higher contrast*/
            var prop = "background-color";
            if ($(this).is('.green')) {
                var v = Math.abs(parseFloat($(this).text().split('%')[0]));
                $(this).css(prop,"rgb(0,"+(Math.round(start+((max-start)*v/max_g))).toString()+",0)");
            } else if ($(this).is('.red')) {
                var v = parseFloat($(this).text().split('%')[0]);
                $(this).css(prop,"rgb("+(Math.round(start+((max-start)*v/max_r))).toString()+",0,0)");
            }
        });

    };
    
    XPT_u.load_graph3 = function ( flag ) {
        XPT_u.create_graph('/dgraph3.png?','cost5','tbi','sp1', flag);
    };

    XPT_u.load_graph4 = function ( flag ) {
        XPT_u.create_graph('/dgraph4.png?','cost5','mdr','sp2', flag);
    };
 
    XPT_u.load_table3 = function () {
        XPT_u.create_table ( 'tbi', 'sp1_table', 'B', 0 );
    };
    
    XPT_u.load_table4 = function () {
        XPT_u.create_table ( 'mdr', 'sp2_table', 'B', 0 );
    };

    XPT_u.load_barbox = function () {
        XPT_u.create_barbox('/bargraph_uncert.png?','bar_box','B');
    };

    XPT_u.load_comp_table = function () {
        XPT_u.create_comp_table( $('#table_box') );
    };

    /* end data load functions */

    XPT_u.unhide_to_start = function () {
        $('#graph_box').show();
        $('#graph1_frame').show();
        $('#linkbar').find('button').each(function () {
            $(this).show();
        });    
    };

    XPT_u.reset_radio = function () {
        $('input:radio[name=sp1_grp]')[0].checked = true;
        /*$('input:radio[name=sp2_grp]')[0].checked = true; Removed until we can figure MDR out */
    };

    /*Event Handling*/

    XPT_u.handle_radio = function ( jqry ) {
        var id = jqry.attr('name').split('_')[0];
        var choice = jqry.val();
        if (id === 'sp1') {
            XPT_u.create_graph('/dgraph3.png?', 'cost5', 'tbi', 'sp1', choice);
        } else if (id === 'sp2') {
            XPT_u.create_graph('/dgraph4.png?', 'cost5', 'mdr', 'sp2', choice);
        }
    };

    $("[name='sp1_grp']").change(function () {
        XPT_u.handle_radio ( $(this) );
    });

    $("[name='sp2_grp']").change(function () {
        XPT_u.handle_radio ( $(this) );
    });

    /* Tooltip Code */

    XPT_u.load_events = function () {
        
        $('[id^="rb"]').click(function () {
            var parse = $(this).attr('id').split('_');
            var type = parse[1];
            var order_by = parse[2];
            var sp_table = (type === 'tbi')?'sp1_table':'sp2_table';
            var radio_selected = $('[name="' + sp_table.split('_')[0] +'_grp"]:checked').val();
            XPT_u.create_table( type, sp_table, radio_selected, order_by );
        });


    };

    /* Initialize */
    
    XPT_u.unhide_to_start();
  //  XPT_u.load_graph3();
   XPT_u.load_table3();

    XPT_u.reset_radio();
    XPT_u.load_name();
    XPT_u.load_m_options();
    XPT_u.load_c_options();
//    XPT_u.load_graph3();
//    XPT_u.load_graph4();
    XPT_u.load_table4();
    XPT_u.load_barbox();
    XPT_u.load_comp_table();

 //   XPT_u.load_events();

});

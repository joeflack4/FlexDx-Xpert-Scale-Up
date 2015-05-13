$(document).ready ( function () {

    (function() { //Hopefully this will fix the IE bug.
        var method;
        var noop = function () {};
        var methods = [
            'assert', 'clear', 'count', 'debug', 'dir', 'dirxml', 'error',
            'exception', 'group', 'groupCollapsed', 'groupEnd', 'info', 'log',
            'markTimeline', 'profile', 'profileEnd', 'table', 'time', 'timeEnd',
            'timeStamp', 'trace', 'warn'
        ];
        var length = methods.length;
        var console = (window.console = window.console || {});

        while (length--) {
            method = methods[length];

        // Only stub undefined methods.
            if (!console[method]) {
                console[method] = noop;
            }
        }
    }());


    function round_to( val, digits ) {
        return Math.round(val * Math.pow(10,digits)) / Math.pow(10,digits);
    }

    var XPT_mdl = {};

    XPT_mdl.wp = $('#web_path').text();
    XPT_mdl.select = parseInt($('#int_select').text(),10);
    XPT_mdl.progress = -1;
    XPT_mdl.json_timer = 5000;
    XPT_mdl.blink_timer = 1000;
    XPT_mdl.graph_toggle = [0,0];
    XPT_mdl.hb = ($("#homebrew_inputs").text() === '') ? false:true;
 
    XPT_mdl.loadInputs = function () {
       $('#in_target_inc').text(XPT_mdl.data.target_inc);
       $('#in_target_hiv').text(XPT_mdl.data.target_hiv);
       $('#in_target_mdr').text(XPT_mdl.data.target_mdr.toString().substring(0,5));
       $('#in_drug1_cost').text(XPT_mdl.data.drug1_cost);
       $('#in_drug2_cost').text(XPT_mdl.data.drug2_cost);
       $('#in_drug3_cost').text(XPT_mdl.data.drug3_cost);
       $('#in_outpt_cost').text(XPT_mdl.data.outpt_cost);
       $('#in_sm_cost').text(XPT_mdl.data.sm_cost);
       $('#in_gxp_cost').text(XPT_mdl.data.gxp_cost);
       $('#in_sdgxp_cost').text(XPT_mdl.data.sdgxp_cost); 
    };

    XPT_mdl.loadData = function ( ) {
        $.ajax ({
            url:XPT_mdl.wp,
            async:false,
            dataType:'json',
            cache:false,
            success: function (data) {
                if (!data.drug1_cost) {
                    //console.log("Data isn't loaded yet, setTimeout()");
                    setTimeout( XPT_mdl.loadData, 500 );
                } else {
                    //console.log("Data Loaded");
                    XPT_mdl.data = data;
                    XPT_mdl.loadInputs();
                }
                //XPT_mdl.select = XPT_mdl.data['int_select'];
            },
            error: function (jq, tx, tr) {
                //console.log('Initial read to JSON file failed',  tx, tr);
                return undefined;
            }
        });
    };
    
    XPT_mdl.loadData();

    (function (mouseOverfile) {
        $.ajax ({
            url:'/static/mouseOver.json',
            async:false,
            dataType:'json',
            cache: false,
            success: function (data) {
                XPT_mdl.mOver = data;
            },
            error: function (jq, tx, tr) {
                return undefined;
                //console.log ('Initial read to MouseOver JSON file failed', tx, tr);
            }
         });
    })();

    $("#left_box").hide();
 
    XPT_mdl.updatePage = function () {
        var var_text = [
                        {'div-cat':'TB Incidence'},
                        {'inc_n':['New:','',' per 100,000']},
                        {'inc_r':['Retreatment:','',' per 100,000']},
                        {'inc_t':['Total:','',' per 100,000']},
                        {'diff-b':['','% ']},
                        {'div-cat':'Drug-Resistance'},
                        {'inc_inh_n':['INH New:','','%']},
                        {'inc_inh_r':['INH Retreatment:','','%']},
                        {'inc_mdr_n':['MDR New:','','%']},
                        {'inc_mdr_r':['MDR Retreatment:','','%']},
                        {'inc_mdr_t':['Total MDR:','',' per 100,000']},
                        {'diff-b':['','% ']},
                        {'div-cat':'Additional Outputs/Indicators'},
                        {'inc_tb_hiv':['TB/HIV:','','%']},
                        {'hiv_prev':['HIV Prevalence:','','%']},
                        {'tb_dur':['TB Duration:','',' years']},
                        {'tb_mort':['TB Mortality:','',' per 100,000']},
                        {'diff-b':['','% ']},
                        {'div-cat':'Costs'},
                        {'cost1':['In Year 1:','$','']},
                        {'diff-b':['','% ']},
                        {'cost5':['In Year 5:','$','']},
                        {'diff-b':['','% ']},
                        {'inter-list':['Year 1','year1']},
                        {'inter-list':['Year 2','year2']},
                        {'inter-list':['Year 3','year3']},
                        {'inter-list':['Year 4','year4']},
                        {'inter-list':['Year 5','year5']},
                        {'inter-list-total':'Total'},
                       ];
        
        var str, insert, prev_key, colour_line, key, direct, interclick_scen, inter_total;
        var total_keys = ['year1','year2','year3','year4','year5'];
        var i, j, k, l, diff, mult, index;
        
        for (j = XPT_mdl.progress; j < XPT_mdl.data.progress; j++) {
            str ="<table class='table table-condensed'>";
            str += "<tr><h4>";
            index = j+1;
            if (XPT_mdl.select != 9 && j != -1) {
                index = XPT_mdl.select;
            }
            if (j != -1) {
                str += (index).toString() + '. ';
            }
            str += XPT_mdl.data[(index).toString()]['name'];
            str += "<tr><th>Description</th><th>Value</th></tr>"
            for (i = 0; i < var_text.length; i ++) {
                interclick_scen = -1;
                key = Object.keys(var_text[i])[0]
                if (key === 'diff-b') {
                    if (j != -1) {
                        prev_key = Object.keys(var_text[i-1])[0];
                        diff = XPT_mdl.data[(index).toString()][prev_key] / XPT_mdl.data['0'][prev_key] * 100.0;
                        diff = 100.0 - diff;
                        if (diff < 0) {
                            colour_line = 'danger';
                            mult = -1;
                            direct = ' increase';
                        } else {
                            colour_line = 'info';
                            mult = 1;
                            direct = ' decrease';
                        }
                        str += "<tr class='" + colour_line + "'><td></td>";
                        str += "<td>";
                        str += (round_to(diff,1) * mult).toString() + var_text[i][key][1];
                        str += direct;
                        str += "</td></tr>";
                    } else {
                        str += "<tr class='active'><td></td><td>Baseline Reference</td></tr>";
                    }
                } else if (key === 'div-cat') {
                    str += "<tr><td colspan='2'><h4>" + var_text[i][key] + "</h4></td></tr>";
                } else if (key === 'inter-list') {
                    if (XPT_mdl.data['homebrew'] === false) { //Not for use in homebrew model
                        if (var_text[i][key][1] === 'year1') {
                            str += "<tr><td colspan='2'><h4>Intermediate Values</h4></td></tr>";
                        }
                        str += "<tr name='interclick_" + (index) + "_" + var_text[i][key][1] + "'>";
                        str += "<td colspan='2'><h5>" + var_text[i][key][0] + "     ";
                        str += "<span class='glyphicon glyphicon-chevron-down'></span></h5></td></tr>";
                        for (k = 0; k < XPT_mdl.data['intermed'].length; k ++) {
                            str += "<tr class='inter_box inter_" + (index) + "_" + var_text[i][key][1] +"'>";
                            str += "<td>" + XPT_mdl.data['intermed'][k] + "</td>";
                            str += "<td>" + XPT_mdl.data[(index).toString()][var_text[i][key][1]][k] + "</td></tr>";
                        }
                    }
                } else if (key === 'inter-list-total') {
                    if (XPT_mdl.data['homebrew'] === false) { //Not for use in homebrew model
                        str += "<tr name='interclick_" + (index) + "_" + var_text[i][key].toLowerCase() + "'>";
                        str += "<td colspan='2'><h5>" + var_text[i][key] + "     ";
                        str += "<span class='glyphicon glyphicon-chevron-down'></span></h5></td></tr>";

                        for (k = 0; k < XPT_mdl.data['intermed'].length; k++) {
                            inter_total = 0;
                            for (l = 0; l < total_keys.length; l++) { //Compute the totals
                                inter_total += XPT_mdl.data[(index).toString()][total_keys[l]][k];
                            }
                            //Totals are computed
                            str += "<tr class='inter_box inter_" + (j+1) + "_" + var_text[i][key].toLowerCase() + "'>";
                            str += "<td>" + XPT_mdl.data['intermed'][k] + "</td>";
                            str += "<td>" + round_to(inter_total,4) + "</td></tr>";
                        }

                        interclick_scen = j+1; //Adds all the events to the dropdowns for this scenario (see below)
                    }

                                    
                } else {
                    str += "<tr><td>" +  var_text[i][key][0] + "</td>";
                    str += "<td>" + var_text[i][key][1] + XPT_mdl.data[(index).toString()][key];
                    str += var_text[i][key][2] + "</td></tr>";
                }
            }
            str += "</table>";
            if (j === -1) {
                $('#left_box div').append(str);
                $("#left_box").show();
            } else {
                if (XPT_mdl.select != 9) {
                    $("#strat" + index.toString() + " div").append(str);
                    $("#strat" + index.toString()).show();
                } else {
                    $('#strat' + (index).toString() + ' div').append(str);
                    $('[name="strat' + (index).toString() + '"]').text(index.toString() + ". " 
                                                                       + XPT_mdl.data[index.toString()].name);
                    $('[name="strat' + (index).toString() + '"]').fadeIn();
                    if (j == 0) {
                        $('#strat1').show();
                    }
                }
            }
            if (interclick_scen != -1) {
                /*Bind the event to the newly added scenario*/
                $("[name^='interclick_" + index + "_']").click( function () {
                    var id = ".inter_" + $(this).attr('name').split('interclick_')[1];
                    $(this).find('span').toggleClass('glyphicon-chevron-down glyphicon-chevron-up');
                    $(id).each( function ( index ) {
                        $(this).slideToggle();
                    });
                });
            }        

        }
        XPT_mdl.progress = XPT_mdl.data.progress;
        XPT_mdl.progress_bar();
    };
    
    XPT_mdl.progress_bar = function () {
        var tot = 10;

        if (XPT_mdl.select != 9) {
            tot = 2;
        } 

        if (XPT_mdl.progress != -1) {
            $("ul#loading li:lt("+(XPT_mdl.progress+1).toString()+")").each(function () {
                $(this).css('background-color', 'rgb(28,72,130)');
            });
        }
        $("#remain").text((tot - XPT_mdl.progress - 1).toString());
       
        if (parseInt($("#remain").text(),10) <= 0) {
            $("#l_console").fadeOut(2000);
        }
    };

    $("[name^='strat']").click( function () {
        var id =$(this).attr('name').split('strat')[1];
        $("[id^='strat']").each ( function (d) {
            $(this).hide();
        });
        $("#strat" + id).show();
    });

    XPT_mdl.getGraphTables = function ( initial ) {
        var i, j, ref_color;
        var t_gph = [];
        var r_gph = [];
        var baseline_x;
        var baseline_y;
        var diff = [ ['TB', 'MDR'], ['tbi','mdr'], ['a','b'], ['inc_t','inc_mdr_t'] ];

        r_gph[0] = XPT_mdl.graph_toggle[0];
        r_gph[1] = XPT_mdl.graph_toggle[1];

        baseline_x = [XPT_mdl.data[r_gph[0].toString()].cost5, XPT_mdl.data[r_gph[1].toString()].cost5];
        baseline_y = [XPT_mdl.data[r_gph[0].toString()].inc_t, XPT_mdl.data[r_gph[1].toString()].inc_mdr_t];
 
        for (i = 0; i < 2; i ++) {
            t_gph[i] = "<tr><th colspan='4'>Percent change in cost and ";
            t_gph[i] += diff[0][i] + " incidence compared to Baseline (Smear) at year 5</th>";
            t_gph[i] += "<th class='center'>Ref. Std.</th></tr>";
            t_gph[i] += "<tr><th>Number</th><th>Name</th><th>Increase in Cost</th><th>Decrease in ";
            t_gph[i] += diff[0][i] + " Incidence</th>";
            t_gph[i] += "<th></th></tr>";

            for (var j=0;j<9;j++) {
                ref_color = '#000000';
                if (r_gph[i] === j) { 
                    ref_color = '#FF0000';
                }
                t_gph[i]+="<tr  style='color:"+ref_color+";'><th class='center'>"+ j.toString()+"</th>";
                t_gph[i]+="<td>" + XPT_mdl.data[j.toString()].name + "</td>";
                t_gph[i] += "<td>"+(round_to((XPT_mdl.data[j.toString()].cost5/baseline_x[i]*100)-100.0,2)).toString()+
                         " %</td><td>"+(round_to(100.0-(XPT_mdl.data[j.toString()][diff[3][i]]/baseline_y[i]*100),2)).toString()+" %</td>";
                          
                t_gph[i] +="<td class='center'><button class='btn btn-sm' id='rb_";
                t_gph[i]+= diff[1][i] + "_"+j+"' type='radio' value='"+j+"'>Ref.</button></td></tr>";
            }

        }

        return [r_gph, t_gph[0], t_gph[1]]

    };

    XPT_mdl.comparisonChart = function ( ) {
        
        var max_g = 0.0;
        var max_r = 0.0;
        var h_a = ['Chg in total inc.',
                   'Chg in MDR inc.',
                   'Chg in TB mort.',
                   'Cost Chg Yr1',
                   'Cost Chg Yr5'];
        var v_names = ['inc_t','inc_mdr_t','tb_mort','cost1','cost5'];
        var t_string = ""

        t_string+="<table class='table disp'>";
        t_string+="<tr><th></th><th>Name";
        for (var j=0;j<5;j++) {
            t_string +="</th><th>"+h_a[j];
        }
        t_string += "</th></tr>"
        for (var j=1;j<9;j++) {
            t_string+="<tr>";
            t_string+="<th class='center'>"+j.toString()+"</th><th>"+XPT_mdl.data[j.toString()].name+"</th>";
            for (var x=0;x<5;x++) {
                var data_value = (XPT_mdl.data[j.toString()][v_names[x]] / XPT_mdl.data['0'][v_names[x]] * 100.0)-100;
                if (data_value < 0) {
                    if (Math.abs(data_value) > max_g) {
                        max_g = Math.abs(data_value);
                    }
                    t_string+="<td class='green'>" + round_to(data_value*-1,1).toString() + "%</td>";  
                } else {
                    if (data_value > max_r) {
                        max_r = data_value;
                    }
                    t_string+="<td class='red'>"+ round_to(data_value,1).toString() +"%</td>";
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

        return {'string':t_string,'max_r':max_r,'max_g':max_g};

       
    };

    XPT_mdl.graphChart = function ( initial ) {
        var gph1, grp2, bar, i, j, diff_list, tmp, cchart_data, tb;
        var t_gph1, t_gph1, r_gph1, r_gph2, gph_info;
        
        gph_info = XPT_mdl.getGraphTables( initial );
        r_gph1 = gph_info[0][0];
        r_gph2 = gph_info[0][1];
        t_gph1 = gph_info[1];
        t_gph2 = gph_info[2];
       
        gph1 = '/dgraph1.png?';
        for (i = 0; i < 9; i ++ ) {
            gph1 += 'c' + i.toString() + '=';
            gph1 += (XPT_mdl.data[i.toString()].cost5 / XPT_mdl.data['0'].cost5 * 100.0) - 100.0 + '&';
            gph1 += 'i' + i.toString() + '=';
            gph1 += 100 - (XPT_mdl.data[i.toString()].inc_t / XPT_mdl.data['0'].inc_t * 100.0) + '&';
        }
        
        gph1 += "ref=" + r_gph1 + "&";
       
        gph1 += "s8name=" + XPT_mdl.data['8'].name;
        /* Insert Table */

      
        if ( !initial ) {
            $("#graph1_frame table").remove();
            $("#graph1_frame").find('img').attr('src',gph1);
        } else {
            $("#sp1").append('<a href="' + gph1 + '" class="btn btn-sm btn-primary d_plot_button" download="inc_cmore5.png">Download Plot</a>');
            $("#graph1_frame > div").append("<img src='" + gph1 + "' alt='Scatterplot1' class='graph' >");
        }
        $("#graph1_frame").append("<table class='table table_hover table_condensed'>" + t_gph1 + "</table>");


        gph2 = '/dgraph2.png?';
        for (i = 0; i < 9; i ++ ) {
            gph2 += 'c' + i.toString() + '=';
            gph2 += (XPT_mdl.data[i.toString()].cost5 / XPT_mdl.data['0'].cost5 * 100.0) - 100.0 + '&';
            gph2 += 'm' + i.toString() + '=';
            gph2 += 100 - (XPT_mdl.data[i.toString()].inc_mdr_t / XPT_mdl.data['0'].inc_mdr_t * 100.0) + '&';
        }
        
        gph2 += "ref=" + r_gph2 + "&";
        gph2 += "s8name=" + XPT_mdl.data['8'].name;

        if ( !initial ) {
            $("#graph2_frame table").remove();
            $("#graph2_frame").find('img').attr('src',gph2);
        } else {
            $("#sp2").append('<a href="' + gph2 + '" class="btn btn-sm btn-primary d_plot_button" download="mdr_cmore5.png">Download Plot</a>');            
            $("#graph2_frame > div").append("<img src='" + gph2 + "' alt='Scatterplot2' class='graph' >");
        }
        $("#graph2_frame").append("<table class='table table_hover table_condensed'>" + t_gph2 + "</table>");
        
        /* Insert Table */
        if ( initial ) {
            diff_list = ['inc_t','inc_mdr_t','tb_mort','cost1','cost5'];

            bar = '/bargraph.png?';
            for (i = 0; i < 5; i ++ ) {
                for (j = 1; j < 9; j ++ ) {
                    bar+='d' + j.toString() + i.toString() + '=';
                    bar+=(XPT_mdl.data['0'][diff_list[i]]/XPT_mdl.data[j.toString()][diff_list[i]]*100.0)-100.0;
                    bar+='&';
                }
            }

            bar += "s8name=" + XPT_mdl.data['8'].name;

            $("#bar_box").append("<img src='" + bar + "' alt='BarGraph' class='graph' >");
            $("#bar_box").prepend("<a href='" + bar + "' class='btn btn-sm btn-primary d_plot_button' download='barchart.png'>Download Plot</a>");

            cchart_data = XPT_mdl.comparisonChart();
            
            tb = $("#table_box");

            tb.append(cchart_data.string);

            tb.find("td").each(function () {
                var start = 50; /*96*/
                var max = 255;  /*230*/ /*Changed to provide higher contrast*/
                var prop = "background-color";
                if ($(this).is('.green')) {
                    var v = Math.abs(parseFloat($(this).text().split('%')[0]));
                    $(this).css(prop,"rgb(0,"+(Math.round(start+((max-start)*v/cchart_data.max_g))).toString()+",0)");
                } else if ($(this).is('.red')) {
                    var v = parseFloat($(this).text().split('%')[0]);
                    $(this).css(prop,"rgb("+(Math.round(start+((max-start)*v/cchart_data.max_r))).toString()+",0,0)");
                }
            });


            XPT_mdl.progress += 1;
            $("[name='strat9']").fadeIn();
            $("[name='strat10']").fadeIn();
            $("[name='strat11']").fadeIn();
            $("[name='data_ex']").fadeIn();


            XPT_mdl.progress_bar();

        }
        /*
        XPT_mdl.load_tooltips();
        */

        $('[id^="rb"]').click(function () {
            var which;
            which = $(this).attr('id').split('rb_')[1].split('_');
            if (which[0] === 'tbi') {

                if (parseInt(which[1],10) !== XPT_mdl.graph_toggle[0]) {
                    XPT_mdl.graph_toggle[0] = parseInt(which[1],10);
                    XPT_mdl.graphChart ( false );
                }
            } else {
                if (parseInt(which[1],10) !== XPT_mdl.graph_toggle[1]) {
                    XPT_mdl.graph_toggle[1] = parseInt(which[1],10); 
                    XPT_mdl.graphChart ( false );
                }
            }
            //console.log(which);
        });
    };

    XPT_mdl.checkProgress = function () {
        $.ajax ({
            url:XPT_mdl.wp,
            dataType:'json',
            cache:false,
            success: function (data) {
                if (data.progress > XPT_mdl.progress) {
                    XPT_mdl.data = data;
                    XPT_mdl.updatePage();
                }
                if (XPT_mdl.select < 9) {
                    if (XPT_mdl.progress != 1) {
                        setTimeout( XPT_mdl.checkProgress, XPT_mdl.json_timer );
                    } 
                } else {
                    if (XPT_mdl.progress != 8) {
                        setTimeout( XPT_mdl.checkProgress, XPT_mdl.json_timer );
                    } else {
                        XPT_mdl.graphChart( true );
                    }
                }
            },
            error: function (jq, tx, tr) {
                return undefined;
                //console.log('Read to JSON file failed', tx, tr);
            }
        });
    };

    XPT_mdl.blink = function () {
        
        var blinker, val;
        blinker = $("ul#loading li:eq(" + (XPT_mdl.progress+1).toString() + ")");
        //Firefox uses spaces, IE doesn't
        if (blinker === undefined) return;

        val = blinker.css('background-color');
        if (val === undefined) return;
        val = val.replace(/\s/g, "");
        if ( val != 'rgb(222,222,222)' ) {
            blinker.css('background-color','rgb(222, 222, 222)');
        } else {
            blinker.css('background-color','rgb(73, 109, 155)');
        } 
        if (XPT_mdl.select < 8) {
            if (XPT_mdl.progress <= 1) {
                setTimeout ( XPT_mdl.blink, XPT_mdl.blink_timer );
            } 
        } else {
            if (XPT_mdl.progress <= 8) { 
                setTimeout ( XPT_mdl.blink, XPT_mdl.blink_timer );     
            } 
        }
    };  
    
    /*Data Export button*/
    $("[name='data_ex']").click ( function () {
        var filename = XPT_mdl.wp.split('/media/')[1];
        console.log(filename);
        location.replace ('/excel/?filename=' + filename );
    });


    setTimeout (XPT_mdl.checkProgress, XPT_mdl.json_timer );
    setTimeout (XPT_mdl.blink, XPT_mdl.blink_timer );
});


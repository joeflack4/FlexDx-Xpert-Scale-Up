$(document).ready(function () {
   
    var XPT_I = {};

    XPT_I.def_input_vals = [['id_t_inc',250], ['id_t_mdr',3.7], ['id_t_hiv',0.83],
                               ['id_t_drug1_cost',500],['id_t_sm_cost',2],
                               ['id_t_gxp_cost',15],['id_t_cx_cost',20],
                               ['id_t_mods_cost',5],['id_t_dst_cost',40],
                               ['id_t_sd_cost',10],['id_t_sdgxp_cost',30],
                               ['id_t_drug2_cost',1000],['id_t_drug3_cost',5000],
                               ['id_t_outpt_cost',10]];
  
    XPT_I.data_interface_map = [['id_t_inc','tb_point'],
                                   ['id_t_mdr','mdr_point'],
                                   ['id_t_hiv','aids_point'],
                                   ['id_t_drug1_cost','drug1_cost'],
                                   ['id_t_outpt_cost','outpt_cost']];
   
    XPT_I.unhide = function () { 
        $("#left_box").show();
        $("#right_box").show(); 
        $("#linkbar li").show();
    };
          
    function round_to( val, digits ) {
        return Math.round(val * Math.pow(10,digits)) / Math.pow(10,digits);
    } 

    function stopRKey(evt) {
        /* Code to prevent the enter key from submitting the form */
        var evt = (evt) ? evt : ((event) ? event : null);
        var node = (evt.target) ? evt.target : ((evt.srcElement) ? evt.srcElement : null);
        if ((evt.keyCode == 13) && (node.type=="text"))  {return false;}
    }

    document.onkeypress = stopRKey; 

    XPT_I.load_mouseOver = function () {
        $.ajax ({
                /*
                  Load the JSON data into a global variable
                  Reading the file only once per page load.
                */
                url: '/static/mouseOver.json',
                async: false,
                dataType: 'json',
                success: function (data) {
                    XPT_I.m_data = data;
                    XPT_I.load_descriptions();
                }
        });
    };

    XPT_I.load_country_input = function ( which ) {
        $.ajax ({
            url:'/static/data_Sept3/' + which + '.json',
            async:true,
            dataType:'json',
            success: function (data) {
                var new_vars = [], i, check;
                var outpt = parseFloat(data['outpt_cost']);
                var drug1 = parseFloat(data['drug1_cost']);
                var sm, gxp;

                for (i = 0; i < 5; i ++ ) {
                    check = parseFloat(data[XPT_I.data_interface_map[i][1]]);
                    if (isNaN(check) || check <= 0) {
                        check = 0.01;
                    }
                    new_vars.push([XPT_I.data_interface_map[i][0], check]);
                }
                new_vars.push(['id_t_drug2_cost',drug1 * 2]);
                new_vars.push(['id_t_drug3_cost',drug1 * 10]);
                if ( (outpt / 5.0) > 2.0 ) {
                    sm = outpt / 5.0;
                } else {
                    sm = 2.0;
                }
                new_vars.push(['id_t_sm_cost',sm]);
                gxp = 15;
                new_vars.push(['id_t_gxp_cost',gxp]);
                new_vars.push(['id_t_cx_cost',sm * 10]);
                new_vars.push(['id_t_mods_cost',sm * 2.5]);
                new_vars.push(['id_t_dst_cost',sm * 20]);
                new_vars.push(['id_t_sd_cost',sm * 5]);
                new_vars.push(['id_t_sdgxp_cost',(4 * sm) + gxp]);

                for (i = 0; i < new_vars.length; i ++) {
                    new_vars[i][1] = round_to (new_vars[i][1],2);
                    $('#'+new_vars[i][0]).val(new_vars[i][1]);
                }         

            },
            error: function () {
                return undefined;
                //console.log('Country Data failed to load');
            }
        });
    };

    XPT_I.populate_dropdown = function (jsonfile) {
            /*Function to populate the dropdown*/
        $.ajax ({
            url:jsonfile,
            dataType:'json',
            success: function (data) {
                sorted_vk = data.sort();
                $.map(sorted_vk, function (data) {
                    $('#country_select').append("<option value='" + data[1] +"'>" + data[0] + "</option>");
                });
                //$('#country_select').append("<option val='" + + "'>" + + "</option>");
            }
        });
    };

    XPT_I.populate_input_vars = function ( which ) {
        if (which == undefined || which === 'NAC') {
            var i;
            for (i = 0; i < XPT_I.def_input_vals.length; i ++) {
                $('#'+XPT_I.def_input_vals[i][0]).val(XPT_I.def_input_vals[i][1]);
                //console.log(XPT_I.def_input_vals[i]);
            }
        } else {
            XPT_I.load_country_input( which );
        }
    };

    XPT_I.loadCountry = function (path) {
        $.ajax ({
            url:path,
            async:false,
            dataType:'json',
            success:function (data) {
                //console.log(data);
                XPT_I.c_data = data;
            },
            error:function() {
                return undefined;
                //console.log('Country data not loaded');
            }
        });
        return XPT_I.c_data;   
    };

    $('#country_select').change(function () {
        var bool = false;
        var selected = $(this).val();
        /*var csrf;*/
        XPT_I.populate_input_vars( selected );
        if (selected === 'NAC') {
            bool = true;
        }
        $(".run_uncertain").each ( function (index) {
            $(this).prop('disabled', bool);
        });
    });

    $('.run_uncertain').click(function () {
        var selected = $('#country_select').val();
        var csrf;
        if (selected != 'NAC') {
            csrf = $("[name='csrfmiddlewaretoken']").val();
            window.location='/model/?country='+selected+'&csrfmiddlewaretoken='+csrf;
        }
    });

    $("#reset_btn").click(function () {
        $('#country_select').val('NAC');
        $(".run_uncertain").each( function (index) {
            $(this).prop('disabled',true);
        });
    });

    $("[name^='t_']").focusout(function() {
        /*
          Sanity checking on the text-field inputs
        */
        var x = $(this);
        var t = $.trim(x.val());
        x.val(t);
        if (t == '' || isNaN(t)) {
            alert("Input value must exist and be a number");
            x.val('');
            setTimeout(function () {
                x.focus();
            }, 1);
        } else {
            switch ($(this).attr('name')) {
                case 't_sm_cost':
                    if (parseFloat($("[name='t_sd_cost']").val(),10) < 
                        parseFloat($(this).val(),10)) {
                        $("[name='t_sd_cost']").val($(this).val());
                    }
                    break;
                case 't_gxp_cost':
                    if (parseFloat($("[name='t_sdgxp_cost']").val(),10) <
                        parseFloat($(this).val(),10)) {
                        $("[name='t_sdgxp_cost']").val($(this).val());
                    }
                    break;
                case 't_cx_cost':
                    if (parseFloat($("[name='t_dst_cost']").val(),10) <
                        parseFloat($(this).val(),10)) {
                        $("[name='t_dst_cost']").val($(this).val());
                    }
                    
            }
        }
    });

    XPT_I.load_descriptions = function () {
        var build_str;
        var keys = Object.keys(XPT_I.m_data);
        build_str = "<div class='list-group'>";
        $.each(keys, function (index, data) {
            if (index < 9) {
                build_str += "<div class='list-group-item'><h5 class='list-group-item-heading'>";
                build_str += XPT_I.m_data[data][0].substring(3);
                build_str += "</h5><p class='list-group-item-text'>";
                build_str += XPT_I.m_data[data][1] + "</p></div>";
            }
        });
        build_str += "</div>";
        //console.log(keys);
        $(".desc_strat_content > div").append(build_str);
    };
    

    $(".desc_strat_title").click (function () {
        $(".desc_strat_content").slideToggle();
        $("#desc_strat_icon").toggleClass("glyphicon-chevron-down glyphicon-chevron-up");
        $("#desc_box").toggleClass("desc_box_visible desc_box_scroll");
    });

    /* Mouseover code to be removed */
/*
    XPT_I.loadSelected = function () {
        
  //        Display the checked strategy in the mouse over <div>
        
        var single_tmp = $("input[name='int_select']:checked").val();
        var a_s_tmp = $("input[name='type_select']:checked").val();
        if (a_s_tmp == '0') { //Single is clicked
            var wh = 'opt' + single_tmp;
            $("#m_over").html("<h3>"+XPT_I.m_data[wh][0]+"</h3>"+XPT_I.m_data[wh][1]);
        } else { // All is clicked
            $("#m_over").html("<h3>"+XPT_I.m_data['opt9'][0]+"</h3>"+XPT_I.m_data['opt9'][1])
        }
    };

    $("[id^='opt']").mouseenter(function() {
        
    //      Display the strategy under the mouse in the mouse over <div>
        
        var wh = $(this).attr('id');
        $("#m_over").html("<h3>"+XPT_I.m_data[wh][0]+"</h3>"+XPT_I.m_data[wh][1]);
    });

    $("[id^='opt']").mouseleave(function () {
        
        //  If the user is no longer hovering over a strategy, default to selected
        
        XPT_I.loadSelected();
    });
    */
    $("[class='single']").click(function () {
        $("#inputs").show();
        //XPT_I.loadSelected();
    });

    $("[class='all']").click(function () {
        $("#inputs").hide();
        //XPT_I.loadSelected();
    });
    /*
    $("[name='int_select']").click(function () {
        
        //  If the user has clicked a strategy, display it
        
        XPT_I.loadSelected();
    }); */

    $("[value='Reset']").click(function () {
        $("#inputs").hide();
    });

    XPT_I.unhide(); 

    XPT_I.load_mouseOver();
    XPT_I.populate_dropdown('/static/data_Sept3/country_data.json');
    /*XPT_I.populate_input_vars();*/
    /*XPT_I.populate_input_vars('AFG');*/
});

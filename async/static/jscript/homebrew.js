$(document).ready(function () {
    
    $('#titlebox').click(function () {
        window.location = '/';
    });
    
      
    var FLEXDX_HB = {};

    FLEXDX_HB.pte_change = [true, true]; //first for TB patients, second for not

    FLEXDX_HB.diag_test_data = { 'SM' : { name:'Smear',   full_name:'Smear', sp_sens: 100.0, sn_sens: 0.0, spec: 99.0,
                                          text: '' },
                                 'GX' : { name:'GXP',     full_name:'GeneXpert', sp_sens: 100.0, sn_sens: 70.0, spec: 99.0,
                                          text: '' },
                                 'CX' : { name:'CXR',     full_name:'Chest X-Ray', sp_sens: 100.0, sn_sens: 95.0, spec: 80.0,
                                          text: 'Used as a screening test only; must confirm with another test.' },
                                 'CU' : { name:'Culture', full_name:'Culture', sp_sens: 100.0, sn_sens: 85.0, spec: 99.0,
                                          text: '' },
                                 'LA' : { name:'LAM',     full_name:'Urine LAM', sp_sens : 50.0, sn_sens: 50.0, spec: 99.0,
                                          text: 'HIV+ only' },
                                 'U1' : { name:'User Test1' },
                                 'U2' : { name:'User Test2' },
                                 'U3' : { name:'User Test3' },
                                 'NO' : { name:'None' },
                                 'LP' : { name:'LPA'  },
                                 'UD' : { name:'User Defined' } };

    FLEXDX_HB.active_sel = undefined;
    FLEXDX_HB.active_menu = undefined;

    FLEXDX_HB.state = { 'hnn' : [ [ 'SM', 'NO' ], ['NO', 'NO' ] ],
                        'hpn' : [ [ 'SM', 'NO' ], ['NO', 'NO' ] ],
                        'hnp' : [ [ 'SM', 'NO' ], ['NO', 'NO' ] ],
                        'hpp' : [ [ 'SM', 'NO' ], ['NO', 'NO' ] ] };

    FLEXDX_HB.current_ud = [ 'User Defined', 'User Test1', 'User Test2', 'User Test3' ];


    /* Allow the button menu to disappear */
    FLEXDX_HB.not_menu_hover = 1;

    $(".test_sel_menu").hover ( function () {
        FLEXDX_HB.not_menu_hover ^= 1;
    });

    $(document).on('mouseup keyup', function( e ){
        if ( FLEXDX_HB.not_menu_hover || e.which==27) 
            if ( FLEXDX_HB.active_menu != undefined )
                FLEXDX_HB.active_menu.hide();
    });
    /* end */

    FLEXDX_HB.load_dtd_content = function () {
        var build_str;
        var keys = Object.keys(FLEXDX_HB.diag_test_data);
        build_str = "<div class='list-group'>";
        $.each(keys, function (index, data) {
            if (data !== 'U1' && data !== 'U2' && data != 'U3' &&
                data !== 'NO' && data !== 'LP' && data != 'UD') {
                build_str += "<div class='list-group-item'><h6 class='list-group-item-heading'>";
                build_str += FLEXDX_HB.diag_test_data[data].full_name;
                if (FLEXDX_HB.diag_test_data[data].name !== FLEXDX_HB.diag_test_data[data].full_name) {
                    build_str += "  <small>(" + FLEXDX_HB.diag_test_data[data].name + ")</small>";
                }
                build_str += "</h6><p class='list-group-item-text'>Sensitivity Sm-: ";
                build_str += FLEXDX_HB.diag_test_data[data].sn_sens.toString();
                build_str += "%, Sensitivity Sm+: ";
                build_str += FLEXDX_HB.diag_test_data[data].sp_sens.toString();
                build_str += "%, Specificity: ";
                build_str += FLEXDX_HB.diag_test_data[data].spec.toString();
                build_str += "%</p>";
                build_str += "<p><em><strong>" + FLEXDX_HB.diag_test_data[data].text + "</strong></em></p></div>";           
            }
        });
        build_str += "</div>";
        $(".dtd_content > div").append(build_str);
    };

    FLEXDX_HB.toggle_ud_glyph = function ( selector ) {
        $(selector).toggleClass('glyphicon-chevron-down glyphicon-chevron-up');
    };

    FLEXDX_HB.toggle_diag_disp = function () {
        $(".diag_content").slideToggle();
    };

    FLEXDX_HB.toggle_dtd_disp = function () {
        $(".dtd_content").slideToggle();
    };

    FLEXDX_HB.toggle_ud_disp = function ( which ) {
        $("#ud_content" + which).slideToggle();
    };
  
    FLEXDX_HB.load_diag_test_data = function ( test ) {
        var ud_elems = ['senssmn', 'senssmp', 'spec', 'cost', 'ltfu' ];
        var which = test.substring(1,2);
        var no_error = true;
        //console.log(which);

        test_data = {};
        $.each ( ud_elems, function (i, val) {
            test_data[val] = parseInt($("#" + val + which).val(),10);
            if ( !test_data[val] && test_data[val] != 0) {
                alert("Missing " + val + " data for User Test " + which);
                no_error = false;
            }
        });
        if (no_error) {
            return test_data;
        }
        return false;
    };

    FLEXDX_HB.load_dst_test_data = function () {
        var ud_dst_elems = ['dst_sens', 'dst_spec', 'dst_cost'];

        var tmp_name;

        test_data = {};
        $.each ( ud_dst_elems, function (i, val) {
            tmp_name = val.split('_')[1];
            test_data[tmp_name] = parseInt ($("#" + val).val(),10);
            if (!test_data[tmp_name] && test_data[tmp_name] != 0) {
                alert("Missing " +tmp_name+" data for DST User Test");
                return false;
            }
        });
        return test_data;


    };
  
    FLEXDX_HB.compute_json = function () {
        var diag_elems = [ ['#ds_hn_n_1', '#ds_hn_n_2', '#ds_hn_n_tb', '#ds_hn_n_not'],
                           ['#ds_hp_n_1', '#ds_hp_n_2', '#ds_hp_n_tb', '#ds_hp_n_not'],
                           ['#ds_hn_p_1', '#ds_hn_p_2', '#ds_hn_p_tb', '#ds_hn_p_not'],
                           ['#ds_hp_p_1', '#ds_hp_p_2', '#ds_hp_p_tb', '#ds_hp_p_not'] ];
        var dst_elems = [ 'drg_n_n', 'drg_p_n', 'drg_n_p', 'drg_p_p' ];
        var ud_elems = [ 'senssmn', 'senssmp', 'spec', 'cost', 'ltfu' ];
        var ud_dst_elems = [ 'dst_sens', 'dst_spec', 'dst_cost' ];
        var json_cat = [ 'hiv-new', 'hiv+new', 'hiv-ret', 'hiv+ret' ];

        var json_out = { "diag" : { "hiv+new" : [], "hiv+ret" : [], "hiv-new": [], "hiv-ret" : [] } };
        var i, key, tmp_name, ch, dst_pre = 'DST_';
        var tests_used = {};
        var test_data, good_data = true;
        var err_msg = '';

        json_out['dst'] = { "hiv+new" : [],  "hiv+ret" : [], "hiv-new": [], "hiv-ret" : [] };
        json_out['tests'] = [];
        
        var norm_tests = {};
        var ud_tests = {};
        var i;
        var temp_int, temp_val;
        var rel_map = {'hnn':['hiv-new', '#ds_hn_n_tb', '#ds_hn_n_not'],
                       'hpn':['hiv+new', '#ds_hp_n_tb', '#ds_hp_n_not'],
                       'hnp':['hiv-ret', '#ds_hn_p_tb', '#ds_hn_p_not'],
                       'hpp':['hiv+ret', '#ds_hp_p_tb', '#ds_hp_p_not'] };

        $.each ( Object.keys(rel_map), function (j, map_val) {
            $.each ( FLEXDX_HB.state[map_val], function (index, value) {
                //$.each ( FLEXDX_HB.state['hnn'][index], function (inner_index, inner_value) {
                if (value[0].substring(0,1) !== 'U') {
                    json_out['diag'][rel_map[map_val][0]].push (FLEXDX_HB.diag_test_data[value[0]].name);
                    if (value[0] !== 'NO')
                        norm_tests[FLEXDX_HB.diag_test_data[value[0]].name] = '';
                } else {
                    temp_int = parseInt(value[0].substring(1,2),10);
                    json_out['diag'][rel_map[map_val][0]].push (FLEXDX_HB.current_ud[temp_int]);
                    temp_val = FLEXDX_HB.load_diag_test_data(value[0])
                    if ( temp_val === false ) {
                        good_data = false;
                    } else {
                        norm_tests[dst_pre + FLEXDX_HB.current_ud[temp_int]] = ''
                        ud_tests[FLEXDX_HB.current_ud[temp_int]] = temp_val;
                    }
                }
                if (value[1] !== 'UD') {
                    json_out['dst'][rel_map[map_val][0]].push (dst_pre + FLEXDX_HB.diag_test_data[value[1]].name);
                    if (value[1] !== 'NO')
                        norm_tests[dst_pre + FLEXDX_HB.diag_test_data[value[1]].name] = '';
                } else {
                    json_out['dst'][rel_map[map_val][0]].push (dst_pre + FLEXDX_HB.current_ud[0]);
                    temp_val = FLEXDX_HB.load_dst_test_data();
                    if ( temp_val === false ) {
                        good_data = false;
                    } else {
                        norm_tests[dst_pre + FLEXDX_HB.current_ud[temp_int]] = ''
                        ud_tests[dst_pre + FLEXDX_HB.current_ud[0]] = temp_val;
                    }
                }
            });
            json_out['diag'][rel_map[map_val][0]].push ( parseInt($(rel_map[map_val][1]).val(),10) / 100.0 );
            json_out['diag'][rel_map[map_val][0]].push ( parseInt($(rel_map[map_val][2]).val(),10) / 100.0 );
            if ( parseInt($(rel_map[map_val][1]).val(),10) < parseInt($(rel_map[map_val][2]).val(),10) ) {
                good_data = false;
                err_msg = '"People with TB symptoms but without true TB" cannot be greater than "People with active TB"';
            }
        });

        json_out['tests'] = Object.keys(norm_tests);
        json_out['ud_tests'] = ud_tests;
                if (err_msg) {
            alert(err_msg);
        }
        if (good_data) {
            var hb_name = $("#homebrew_name");
            if (!hb_name.val()) {
                hb_name.val("User Defined");
            }  
            json_out['ud_strat_name'] = hb_name.val();
            $("#json_output > p").text(JSON.stringify(json_out, null, '\t'));
            //console.log(JSON.stringify(json_out, null, '\t'));
            $("#json_object").val(JSON.stringify(json_out));
        }
        //console.log(good_data);
        return good_data;
    };
    
    /* Managing the title clicks for the user defined strats (Show/Hide) */
    $(".testbox_title").click ( function () {
        var which_ud = $(this).attr('id').split('title')[1];
        FLEXDX_HB.toggle_ud_disp( which_ud );
        FLEXDX_HB.toggle_ud_glyph("#test_icon" + which_ud);
    });

    /* Managing the Diagnostic/Treatment/Description box clicks  (Show/Hide) */
    $("#diag_title").click ( function () {
        FLEXDX_HB.toggle_diag_disp();
        FLEXDX_HB.toggle_ud_glyph("#diag_icon");
    });

    $("#dtd_title").click ( function () {
        FLEXDX_HB.toggle_dtd_disp();
        FLEXDX_HB.toggle_ud_glyph('#dtd_icon');
    });


    $("[name='pte_input']").change ( function () {
        var id = $(this).attr('id');
        var hiv = id.split('_')[1];
        var treat = id.split('_')[2];
        var tb_index = (id.split('_').pop() === 'not') ? 1 : 0;
        var tb_index_val = id.split('_').pop();

        var val;
        if (hiv === 'hn' && treat === 'n') {
            if (FLEXDX_HB.pte_change[tb_index]) {
                
                val = $(this).val();
                $("[name='pte_input']").each ( function (index, element) {
                    if (tb_index_val === $(element).attr('id').split('_').pop()) {
                        $(element).val(val);
                    }
                });

                FLEXDX_HB.pte_change[tb_index] = false;
            }
        } else {
            FLEXDX_HB.pte_change[0] = false;
            FLEXDX_HB.pte_change[1] = false;
        }

    });

    /* Save the User Defined Scenario */
    $("[id^='test_save']").click (function () {
        var which_ud = $(this).attr('id').split('save')[1];
        var test_name = $("#test_title" + which_ud).val(); //User Defined Entry
        var ud_name_map = { '1' : [ 'm1a_ud1_name', 'm1b_ud1_name', 'm2a_ud1_name', 'm2b_ud1_name'],
                            '2' : [ 'm1a_ud2_name', 'm1b_ud2_name', 'm2a_ud2_name', 'm2b_ud2_name'],
                            '3' : [ 'm1a_ud3_name', 'm1b_ud3_name', 'm2a_ud3_name', 'm2b_ud3_name'],
                            '4' : ['mn1a_sm_ud','mn1a_gxp_ud','mn1a_cul_ud','mn1a_lam_ud',
                                   'mn1a_ud1_ud','mn1a_ud2_ud','mn1a_ud3_ud',
                                   'mn1b_sm_ud','mn1b_gxp_ud','mn1b_cul_ud',
                                   'mn1b_ud1_ud','mn1b_ud2_ud','mn1b_ud3_ud',
                                   'mn2a_sm_ud','mn2a_gxp_ud','mn2a_cul_ud','mn2a_lam_ud',
                                   'mn2a_ud1_ud','mn2a_ud2_ud','mn2a_ud3_ud',
                                   'mn2b_sm_ud','mn2b_gxp_ud','mn2b_cul_ud',
                                   'mn2b_ud1_ud','mn2b_ud2_ud','mn2b_ud3_ud'] };
        if ( test_name !== "") {
            /* Add it to the title bar */
            $("#test_name" + which_ud).text(test_name); //Entry in titlebar

            if (which_ud !== '4' ) {
                FLEXDX_HB.current_ud [parseInt(which_ud, 10)] = test_name;
            } else {
                FLEXDX_HB.current_ud [ 0 ] = test_name;
            }

            FLEXDX_HB.reload_selectors();
            $.each (ud_name_map[which_ud], function (index, value) {
                $("#"+value).text(test_name);
            });

        }

        //For now assume the test is completed with this click
        $("#ud_title" + which_ud).removeClass('test_empty').addClass('test_completed');

    });

    $("[id^='test_clear']").click (function () {
        var which_ud = $(this).attr('id').split('clear')[1];
        var ud_name_map = { '1' : [ 'm1a_ud1_name', 'm1b_ud1_name', 'm2a_ud1_name', 'm2b_ud1_name'],
                            '2' : [ 'm1a_ud2_name', 'm1b_ud2_name', 'm2a_ud2_name', 'm2b_ud2_name'],
                            '3' : [ 'm1a_ud3_name', 'm1b_ud3_name', 'm2a_ud3_name', 'm2b_ud3_name'],
                            '4' : ['mn1a_sm_ud','mn1a_gxp_ud','mn1a_cul_ud','mn1a_lam_ud',
                                   'mn1a_ud1_ud','mn1a_ud2_ud','mn1a_ud3_ud',
                                   'mn1b_sm_ud','mn1b_gxp_ud','mn1b_cul_ud',
                                   'mn1b_ud1_ud','mn1b_ud2_ud','mn1b_ud3_ud',
                                   'mn2a_sm_ud','mn2a_gxp_ud','mn2a_cul_ud','mn2a_lam_ud',
                                   'mn2a_ud1_ud','mn2a_ud2_ud','mn2a_ud3_ud',
                                   'mn2b_sm_ud','mn2b_gxp_ud','mn2b_cul_ud',
                                   'mn2b_ud1_ud','mn2b_ud2_ud','mn2b_ud3_ud'] };

        $('#test_name' + which_ud).text(''); //Entry in titlebar
        $('#test_title' + which_ud).val(''); //User defined Entry

        if (which_ud !== '4') {
            FLEXDX_HB.current_ud [parseInt(which_ud, 10) ] = FLEXDX_HB.diag_test_data['U' + which_ud].name;
            $.each (ud_name_map[which_ud], function (index, value) {
                $("#"+value).text(FLEXDX_HB.diag_test_data['U' + which_ud].name);
            });
        } else {
            FLEXDX_HB.current_ud[0] = FLEXDX_HB.diag_test_data['UD'].name;
            $.each (ud_name_map[which_ud], function (index, value) {
                $("#"+value).text(FLEXDX_HB.diag_test_data['UD'].name);
            });
        }

        FLEXDX_HB.reload_selectors();
        $("#ud_title" + which_ud).removeClass('test_completed').addClass('test_empty');

    });

    $("#run_model").click( function () {
        var good_data = FLEXDX_HB.compute_json();
        var json_data;
        var good_model_data = true;

        if (good_data) {
            json_data = $.parseJSON($("#json_object").val());
            json_data['model_inputs'] = {};
            $("#epi_model_inputs input").each (function (index) {
                if ( !$.isNumeric ( $(this).val() ) ) { 
                    good_model_data = false;
                    alert ("The model input value " + $(this).attr('id').split('id_t_')[1] + " is not a proper number");
                } else {
                    json_data['model_inputs'][ $(this).attr('id').split('id_t_')[1] ] = $(this).val();
                }
            });
            if ( good_model_data ) {
                $("#json_object").val(JSON.stringify(json_data));
                $("#homebrew_submit").submit();
            }
        }

    });

    $("#json_only").click ( function () {
        FLEXDX_HB.compute_json();
    });

    $("#hb_input_button").click ( function () {
        $("#homebrew_input_box").toggle();
        if ( $("#homebrew_input_box").is(":visible") ) {
            window.location.href = '#diag_title';
        }
    });

    /* Select Menu1 Init */
    $( "#menu1a" ).menu({
        items: ".menu_item"
    });
    $(" #menu2a" ).menu({
        items: ".menu_item"
    });
    $( "#menu1b" ).menu({
        items: ".menu_item"
    });
    $(" #menu2b" ).menu({
        items: ".menu_item"
    });



    $( "#menu1a" ).hide();
    $( "#menu2a" ).hide();
    $( "#menu1b" ).hide();
    $( "#menu2b" ).hide();


    $(".test_selector").click ( function () {

        var pos = $(this).position();
        var which_num = $(this).attr('id').substring(3,4);
        var which_let = ($(this).attr('id').substring(1,2) === 'n')? 'b' : 'a';
        var menustring = "#menu";
        var height = parseInt($(this).css('height'),10);

        $("[id^='menu']").each ( function () {
            $(this).hide();
        });
       
        FLEXDX_HB.active_sel = $(this);
        
        menustring = menustring + which_num + which_let;
        FLEXDX_HB.active_menu = $(menustring);

        $(menustring).css('top',(pos.top + height).toString() + "px");
        $(menustring).css('left',(pos.left + 15).toString() + "px");
        $(menustring).show();
    });


    FLEXDX_HB.reload_selectors = function () {
 
         var test_name, dst_name;
         var map = Object.keys(FLEXDX_HB.state);
         $.each(map, function (index, value) {
            $.each(FLEXDX_HB.state[value], function (inner_index, inner_value) {
                // Test 
                if (inner_value[0].substring(0,1) === 'U') {
                    //Pull the user defined name
                    test_name = FLEXDX_HB.current_ud [ parseInt(inner_value[0].substring(1,2)) ];
                } else {
                    test_name = FLEXDX_HB.diag_test_data[inner_value[0]].name;
                }
                
                $("#" + value + (inner_index + 1) + '_test').text (test_name);
                
                // DST
                if (inner_value[1].substring(0,1) === 'U' ) {
                    //Pull the user defined name
                    dst_name = FLEXDX_HB.current_ud [ 0 ];
                } else {
                    dst_name = FLEXDX_HB.diag_test_data[inner_value[1]].name;
                }
                
                $("#" + value + (inner_index + 1) + '_dst').text (dst_name);
                
            });
         });
    };

    $('[id^="mn"]').click ( function () {
        var test = $(this).attr('id').split('_')[1];
        var dst = $(this).attr('id').split('_')[2];
        var which_num = $(this).attr('id').substring(2,3);
        var which_let = $(this).attr('id').substring(3,4);
        
        var id = FLEXDX_HB.active_sel.attr('id').split('_')[0];
        var wcat = id.substring(0,3); 
        var wtest = parseInt(id.substring(3,4), 10) - 1;

        var name_map = {'sm':'SM', 'none':'NO', 'gxp':'GX', 'lpa':'LP','cul':'CU','cxr':'CX','lam':'LA',
                        'ud':'UD', 'ud1':'U1', 'ud2':'U2', 'ud3':'U3' };

        console.log(wcat, wtest);

        FLEXDX_HB.state[wcat][wtest][0] = name_map[test];
        FLEXDX_HB.state[wcat][wtest][1] = name_map[dst];

        console.log ( FLEXDX_HB.state );

        FLEXDX_HB.reload_selectors();
        
        $( "#menu" + which_num + which_let ).hide(); 
    });

     /*Initialization*/
    FLEXDX_HB.load_dtd_content();

});

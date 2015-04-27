$(document).ready( function () {

    $('#titlebox').click(function () {
        window.location = '/index.html';
    });

    $("[name^='strat']").click(function () { 
        /*When tabs are clicked display: run_all*/
        var wh = $(this).attr('name');
        switch (wh) {
            case 'strat9': /* Strat9 - Scatter Cost/Incidence graph_box > graph1_frame */
                $("#left_box").hide();
                $("[id^='strat']").hide();
                $("#summary_box").hide();
                $("#graph_box").show();
                $("#graph2_frame").hide();
                $("#graph1_frame").show();
                break;
            case 'strat10': /* Strat10 - Scatter Cost/MDR graph_box > graph2_frame */
                $("#left_box").hide();
                $("#summary_box").hide();
                $("[id^='strat']").hide();
                $("#graph_box").show();
                $("#graph1_frame").hide();
                $("#graph2_frame").show();
                break;
            case 'strat11': /* Strat11 - Summary box (Table and Barchart) */
                $("#left_box").hide();
                $('#graph_box').hide();
                $('#bar_box').show();
                $('#table_box').show();
                $("[id^='strat']").hide();
                $("#summary_box").show();
                break;
            default: 
                wh = "#" + wh;
                $("[id^='strat']").hide();
                $("#graph_box").hide();
                $("#bar_box").hide();
                $("#table_box").hide();
                $("#left_box").show();
                $(wh).show();
        }
    }); 

});

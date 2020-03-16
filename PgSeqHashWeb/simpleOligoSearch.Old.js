
var bgCnt = 0;
var sosThings = new Object;

/* Lots of ExtJS stuff is being guided by code written by Mark Russo,
   in some cases taken in whole cloth from it

   http://kraken.pri.bms.com/biohtml/russom/dev/aso/public_html/extjsgrid2.html

   ExtJS / Sencha docs:
   Store:    http://docs.sencha.com/extjs/4.2.1/#!/api/Ext.data.Store
   MemProxy: http://docs.sencha.com/extjs/4.2.1/#!/api/Ext.data.proxy.Memory
   
*/

Ext.QuickTips.init();

var svcCb = function (rsp, atObj) {
    var data = JSON.parse( rsp.responseText );
    if (!data) return;
    
    var st = data.svctype;
    if (!st) {
        atObj.dump_data( data );
        return;
        // Default AjaxTask handling:
        return atObj.set_content( data );
    }
    var store;
    if (st == 'summary_grid') {
        store = summary_store();
    } else if (st == 'detailed_grid') {
        store = detailed_store();
    } else if (st == 'overview_grid') {
        store = overview_store();
    } else {
        atObj.error_msg("Unknown service type '"+st+"'");
    }
    if (store) {
        add_to_store(data, store, atObj);
        remove_task( atObj );
        return;
    } else {
        return atObj.error_msg("Not sure how to handle output");
    }
}

function strand_renderer (val, metaData) {
    var cls;
    if (val == null) {
        cls = 'nostrand';
    } else if (val == 1) {
        cls = 'strF';
    } else if (val == -1) {
        cls = 'strR';
    } else {
        cls = 'strUnk';
    }
    if (cls) metaData.tdCls = cls;
    return val;
}

function mismatch_renderer (val, metaData) {
    var cls;
    if (val == null) {
        cls = 'mmunk';
    } else if (val == 0) {
        cls = 'mm0';
    } else if (val == 1) {
        cls = 'mm1';
    } else if (val == 2) {
        cls = 'mm2';
    } else {
        cls = 'mm3';
    }
    if (cls) metaData.tdCls = cls;
    return val;
}

function ambig_renderer (val, metaData) {
    var cls;
    if (val == null || val == 0) {
        val = '';
    } else {
        metaData.tdCls = cls = 'mmamb';
    }
    return val;
}

function symbol_list_renderer  (val, metaData) {
    var to = typeof(val);
    if (to == 'undefined') return val;
    if (to != 'object') return val;
    var list = new Array;
    var lnum = val.length;
    for (var i = 0; i < lnum; i++) {
        var gdat  = val[i];
        var acc   = gdat[0];
        var title = gdat[2] ? " title='"+quoteattr(gdat[2])+"'" : "";
        var label = gdat[1];
        var cls   = 'sym';
        if (!label) {
            label = acc;
            cls   = 'acc';
        }
        list.push( "<a lid='"+acc+"' onclick='symPop(this)' class='"+cls+"'"+title+">"+label+"</a>");
    }
    return list.join(', ');
}

// http://www.sencha.com/forum/showthread.php?235912-window-issue&p=868167#post868167
// Do not need this - seems to be for a window in a window.
Ext.define('locusPopup',{
    extend:'Ext.window.Window',
    title:'Locus Details',
    height:100,
    width:100,
    modal:true
});

function symPop (anc) {
    if (!anc) return false;
    var lid = anc.getAttribute('lid');
    if (!lid) return false;
    var title = anc.innerHTML;
    if (anc.title) title = title + " : " + anc.title;
    Ext.create('Ext.window.Window', {
        title: title,
        height: 500,
        width: 600,    
        items:[{
            xtype : 'box',
            html:'Really cool things will be here. They will amaze your current friends, win you buckets of new ones, count as your Fifth Powerful Conversation, and guarantee an Exceeding rating, <i>A++++, would hire again</i>.',
        }]
    }).show();
    return true;
}

function exintmm_renderer (val, metaData) {
    var col;
    val = parseInt(val);
    if (val == null || val == 0) {
        col = '#ddd';
        val = null;
    } else if (val == 1) {
        col = '#0f0';
    } else if (val <= 5) {
        col = '#cf0';
    } else if (val <= 10) {
        col = '#ff0';
    } else if (val <= 50) {
        col = '#f90';
    } else if (val <= 100) {
        col = '#f39';
    } else if (val <= 500) {
        // #fff did not work for some reason...
        col = 'red';
    } else {
        col = 'red';
    }
    if (col) metaData.tdAttr = 'bgcolor='+ col;
    return val;
}

function exonic_renderer (val, metaData) {
    if (val == null) {
    } else {
        var all = val.split(/[, ]/);
        metaData.tdCls = all[0].toLowerCase();
    }
    return val;
}

function note_renderer (val, metaData) {
    if (val) {
        // metaData.tdAttr = 'data-qtip="' + 'Hello' + '"';
        // metaData.attr = 'ext:qtip="' + "hello" + '"';
        metaData.tdAttr = 'data-qtip="' + val + '"';
        metaData.tdCls = 'alert';
        return "!";
    }
    return "";
}

function rnafoot_renderer (val, metaData) {
    if (val) {
        var all = val.split(/[\s\,]+/);
        val = all.join('<br />');
    }
    return val;
}


function remove_task (atObj) {
    var targ = atObj.target();
    if (!targ) return;
    var par = targ.parentNode;
    if (!par) return;
    par.removeChild(targ);
    var kids = par.getElementsByTagName('div');
    if (kids.length) return;
    var gpar = par.parentNode;
    if (gpar) gpar.removeChild( par );
}

function launch_service (params, id) {
    if (!params.bgnd) {
        params.bgnd = --bgCnt;
    }
    new AjaxTask( 'simpleOligoSearch.pl', params, id, svcCb );
   
    return;
}



function add_to_store (data, store, atObj) {
    if (!data || !store) return;
    var rows = data.rows;
    if (rows) {
        return store.loadData( rows, true );
    } else {
        return atObj.error_msg("Failed to find data in returned information");
    }
}

function summary_model () {
    var mName = 'SummaryModel';
    var obj = sosThings.summaryModel;
    if (obj) return mName;
    
    obj = sosThings.summaryModel = Ext.define(mName, {
        extend: 'Ext.data.Model',
        fields: [
            {name: 'oligoid',   type: 'string'},
            {name: 'oligoseq',  type: 'string'},
            {name: 'build',     type: 'string'},
            {name: 'mm0',       type: 'int'},
            {name: 'mm1',       type: 'int'},
            {name: 'mm2',       type: 'int'},
            {name: 'database',  type: 'string'}
        ]
    });
    return obj;
}

function summary_store ( gridTarg ) {
    var obj = sosThings.summaryStore;
    if (obj) return obj;
    var model = summary_model();
    obj = sosThings.summaryStore = Ext.create('Ext.data.Store', {
        model: model.getName(),
        proxy: {
            type: 'memory',
            reader: {
                type: 'json',
                root: 'users'
            }
        },
        autoLoad: false
    });
    // Create a single grid at this time
    summary_grid( gridTarg );
    return obj;
}

function summary_grid (targ) {
    var obj = sosThings.summaryTable;
    if (obj) return obj;
    var store = summary_store();
    if (!targ) targ = 'extout'; // Ext.getBody();
    // Taken from Mark Russo's example
    var mmWid = 40;
    obj = sosThings.summaryTable = new Ext.grid.Panel({
        title: 'Summary ASO Results',
        width: 900,                 // Dimensions
        height: 300,
        renderTo: targ,
        collapsible: true,          // collapes using double arrows in header.
        store: store,               // The store object that manages data.
        multiSelect: true,          // Allows the user to select multiple rows
        plugins: 'gridfilters',     // adds a simple and flexible filter set
        
        // Defines the mask to display while data are being loaded.
        // If false, will show no mask. If true, will display the
        // value of loadingText
        viewConfig: {
            loadMask: true, //false
            loadingText: 'Downloading results...'
        },

        // Remove records from the selection when they are removed
        // from the store.  Important: When using paging or a
        // Ext.data.BufferedStore, records which are cached in the
        // Store's data collection may be removed from the Store when
        // pages change, or when rows are scrolled out of view. For
        // this reason pruneRemoved should be set to false when using
        // a buffered Store.
        selModel: {
            pruneRemoved: false
        },
        
        // All column definitions
        columns:[{
            xtype: 'rownumberer',
            width: 50,
            sortable: false
        },{
            text: "OligoID",
            dataIndex: 'oligoid',
            width: 100,
            sortable: true,
            cellWrap: false,
            filter : { type : 'list' }
        },{
            text: "OligoSeq",
            dataIndex: 'oligoseq',
            align: 'left',
            width: 150,
            sortable: true,
            filter: { type: 'string' }
        },{
            text: "Build",
            dataIndex: 'build',
            align: 'center',
            width: 100,
            sortable: true,
            filter: { type: 'list' }
        },{
            text: "MM0",
            dataIndex: 'mm0',
            align: 'center',
            width: mmWid,
            sortable: true,
            filter: { type: 'numeric' },
            renderer: exintmm_renderer
        },{
            text: "MM1",
            dataIndex: 'mm1',
            align: 'center',
            width: mmWid + 5,
            sortable: true,
            filter: { type: 'numeric' },
            renderer: exintmm_renderer
        },{
            text: "MM2",
            dataIndex: 'mm2',
            align: 'center',
            width: mmWid + 10,
            sortable: true,
            filter: { type: 'numeric' },
            renderer: exintmm_renderer
        },{
            text: "Database Path",
            dataIndex: 'database',
            align: 'left',
            width: 300,
            sortable: true,
            filter: { type: 'string' }
        }]
    });
    return obj;
}

function detailed_model () {
    var mName = 'DetailedModel';
    var obj = sosThings.detailedModel;
    if (obj) return obj;
    
    obj = sosThings.detailedModel = Ext.define(mName, {
        extend: 'Ext.data.Model',
        fields: [
            {name: 'oligoid',   type: 'string'},
            {name: 'sym',       type: 'string'},
            {name: 'build',     type: 'string'},
            {name: 'mm',        type: 'int'},
            {name: 'am',        type: 'int'},
            {name: 'str',       type: 'int'},
            {name: 'exon',      type: 'string'},
            {name: 'alignment', type: 'string'},
            {name: 'geneid',    type: 'string'},
            {name: 'genedesc',  type: 'string'},
            {name: 'sbjid',     type: 'string'},
            {name: 'sbjpos',    type: 'int'},
            {name: 'footprint', type: 'string'},
            {name: 'rnafoot',   type: 'string'},
            {name: 'exonnum',   type: 'string'}, // NOT INTEGER!
            {name: 'notes',     type: 'string'}
        ]
    });
    return obj;
}

function detailed_store ( gridTarg ) {
    var obj = sosThings.detailedStore;
    if (obj) return obj;
    var model = detailed_model();
    obj = sosThings.detailedStore = Ext.create('Ext.data.Store', {
        model: model.getName(),
        proxy: {
            type: 'memory',
            reader: {
                type: 'json',
                root: 'users'
            }
        },
        autoLoad: false
    });
    // Create a single grid at this time
    detailed_grid( gridTarg );
    return obj;
}

function detailed_grid (targ) {
    var obj = sosThings.detailedTable;
    if (obj) return obj;
    var store = detailed_store();
    if (!targ) targ = 'extout';
    // Taken from Mark Russo's example
    obj = sosThings.detailedTable = new Ext.grid.Panel({
        title: 'Detailed ASO Results',
        width: 900,                 // Dimensions
        height: 700,
        renderTo: targ,
        collapsible: true,          // collapes using double arrows in header.
        store: store,               // The store object that manages data.
        multiSelect: true,          // Allows the user to select multiple rows
        plugins: 'gridfilters',     // adds a simple and flexible filter set
        
        // Defines the mask to display while data are being loaded.
        // If false, will show no mask. If true, will display the
        // value of loadingText
        viewConfig: {
            loadMask: false, // true, //false
            loadingText: 'Downloading results...'
        },

        // Remove records from the selection when they are removed
        // from the store.  Important: When using paging or a
        // Ext.data.BufferedStore, records which are cached in the
        // Store's data collection may be removed from the Store when
        // pages change, or when rows are scrolled out of view. For
        // this reason pruneRemoved should be set to false when using
        // a buffered Store.
        selModel: {
            pruneRemoved: false
        },
        
        // All column definitions
        columns:[{
            xtype: 'rownumberer',
            width: 50,
            sortable: false
        },{
            text: "OligoID",
            dataIndex: 'oligoid',
            width: 100,
            sortable: true,
            cellWrap: false,
            filter : { type : 'list' }
        },{
            text: "Symbol",
            dataIndex: 'sym',
            width: 100,
            sortable: true,
            filter : { type : 'string' }
        },{
            text: "Build",
            dataIndex: 'build',
            align: 'center',
            width: 100,
            sortable: true,
            filter: { type: 'list' }
        },{
            text: "MM",
            dataIndex: 'mm',
            align: 'center',
            width: 40,
            sortable: true,
            filter: { type: 'list' },
            renderer: mismatch_renderer
        },{
            text: "Str",
            dataIndex: 'str',
            align: 'center',
            width: 40,
            sortable: true,
            filter: { type: 'list' },
            renderer: strand_renderer
        },{
            text: "Exonic",
            dataIndex: 'exon',
            align: 'center',
            width: 100,
            sortable: true,
            filter: { type: 'list' },
            renderer: exonic_renderer
        },{
            text: "!!",
            dataIndex: 'notes',
            align: 'center',
            width: 30,
            sortable: true,
            filter: { type: 'string' },
            getTip: function () { return "hello"; },
            renderer: note_renderer
        },{
            text: "Alignment",
            dataIndex: 'alignment',
            align: 'left',
            width: 200,
            tdCls: 'seq',
            sortable: false
        },{
            text: "Gene Description",
            dataIndex: 'genedesc',
            align: 'left',
            width: 300,
            sortable: true,
            filter: { type: 'string' }
        },{
            text: "GeneId",
            dataIndex: 'geneid',
            align: 'center',
            width: 100,
            sortable: true,
            hidden: true
        },{
            text: "Genomic Footprint",
            dataIndex: 'footprint',
            align: 'left',
            width: 200,
            sortable: true,
            hidden: true,
            filter: { type: 'string' }
        },{
            text: "AM",
            dataIndex: 'am',
            align: 'center',
            width: 40,
            sortable: true,
            hidden: true,
            filter: { type: 'list' },
            renderer: ambig_renderer
        },{
            text: "RNA Footprint",
            dataIndex: 'rnafoot',
            align: 'left',
            width: 200,
            sortable: true,
            hidden: true,
            filter: { type: 'string' },
            renderer: rnafoot_renderer
        },{
            text: "Exon#",
            dataIndex: 'exonnum',
            align: 'center',
            width: 90,
            sortable: true,
            hidden: true,
           filter: { type: 'string' }
        },{
            text: "Subject ID",
            dataIndex: 'sbjid',
            align: 'left',
            width: 200,
            sortable: true,
            hidden: true
        },{
            text: "Subject Pos",
            dataIndex: 'sbjpos',
            align: 'left',
            width: 80,
            sortable: true,
            hidden: true
        }]
    });
    return obj;
}

function overview_model () {
    var mName = 'OverviewModel';
    var obj = sosThings.overviewModel;
    if (obj) return obj;
    
    obj = sosThings.overviewModel = Ext.define(mName, {
        extend: 'Ext.data.Model',
        fields: [
            {name: 'oligoid',   type: 'string'},
            {name: 'oligoseq',  type: 'string'},
            {name: 'build',     type: 'string'},
            
            {name: 'ex0',       type: 'int'},
            {name: 'int0',      type: 'int'},
            {name: 'ex1',       type: 'int'},
            {name: 'int1',      type: 'int'},
            {name: 'ex2',       type: 'int'},
            {name: 'int2',      type: 'int'},
            
            {name: 'ex0sym',    type: 'auto'},
            {name: 'int0sym',   type: 'auto'},
            {name: 'ex1sym',    type: 'auto'},
            {name: 'int1sym',   type: 'auto'},
            {name: 'ex1sym',    type: 'auto'},
            {name: 'int1sym',   type: 'auto'}
        ]
    });
    return obj;
}

function overview_store ( gridTarg ) {
    var obj = sosThings.overviewStore;
    if (obj) return obj;
    var model = overview_model();
    obj = sosThings.overviewStore = Ext.create('Ext.data.Store', {
        model: model.getName(),
        proxy: {
            type: 'memory',
            reader: {
                type: 'json',
                root: 'users'
            }
        },
        autoLoad: false
    });
    // Create a single grid at this time
    overview_grid( gridTarg );
    return obj;
}

function overview_grid (targ) {
    var obj = sosThings.overviewTable;
    if (obj) return obj;
    var store = overview_store();
    if (!targ) targ = 'extout';
    var mmWid = 20;
    var symWid = 150;
    // Taken from Mark Russo's example
    obj = sosThings.overviewTable = new Ext.grid.Panel({
        title: 'Overview ASO Results',
        width: 900,                 // Dimensions
        height: 400,
        renderTo: targ,
        collapsible: true,          // collapes using double arrows in header.
        store: store,               // The store object that manages data.
        multiSelect: true,          // Allows the user to select multiple rows
        plugins: 'gridfilters',     // adds a simple and flexible filter set
        
        // Defines the mask to display while data are being loaded.
        // If false, will show no mask. If true, will display the
        // value of loadingText
        viewConfig: {
            loadMask: false, // true, //false
            loadingText: 'Downloading results...'
        },

        // Remove records from the selection when they are removed
        // from the store.  Important: When using paging or a
        // Ext.data.BufferedStore, records which are cached in the
        // Store's data collection may be removed from the Store when
        // pages change, or when rows are scrolled out of view. For
        // this reason pruneRemoved should be set to false when using
        // a buffered Store.
        selModel: {
            pruneRemoved: false
        },
        
        // All column definitions
        columns:[{
            xtype: 'rownumberer',
            width: 50,
            sortable: false
        },{
            text: "OligoID",
            dataIndex: 'oligoid',
            width: 100,
            sortable: true,
            cellWrap: false,
            filter : { type : 'list' }
        },{
            text: "Build",
            dataIndex: 'build',
            align: 'center',
            width: 100,
            sortable: true,
            filter: { type: 'list' }
        },{
            
            text: "Ex0",
            dataIndex: 'ex0',
            align: 'center',
            width: mmWid,
            sortable: true,
            filter: { type: 'number' },
            renderer: exintmm_renderer
        },{
            text: "Int0",
            dataIndex: 'int0',
            align: 'center',
            width: mmWid,
            sortable: true,
            filter: { type: 'number' },
            renderer: exintmm_renderer
        
        },{
            text: "Ex1",
            dataIndex: 'ex1',
            align: 'center',
            width: mmWid,
            sortable: true,
            filter: { type: 'number' },
            renderer: exintmm_renderer
        },{
            text: "Int1",
            dataIndex: 'int1',
            align: 'center',
            width: mmWid + 10,
            sortable: true,
            filter: { type: 'number' },
            renderer: exintmm_renderer
        },{
            text: "Ex2",
            dataIndex: 'ex2',
            align: 'center',
            width: mmWid + 15,
            sortable: true,
            filter: { type: 'number' },
            renderer: exintmm_renderer
        },{
            text: "Int2",
            dataIndex: 'int2',
            align: 'center',
            width: mmWid + 20,
            sortable: true,
            filter: { type: 'number' },
            renderer: exintmm_renderer
        },{

            
            text: "Exon 0 Genes",
            dataIndex: 'ex0sym',
            align: 'left',
            width: symWid,
            sortable: true,
            filter: { type: 'string' },
            renderer: symbol_list_renderer
        },{
            text: "Intron 0 Genes",
            dataIndex: 'int0sym',
            align: 'left',
            width: symWid,
            sortable: true,
            filter: { type: 'string' },
            renderer: symbol_list_renderer
        },{
            text: "Exon 1 Genes",
            dataIndex: 'ex1sym',
            align: 'left',
            width: symWid,
            sortable: true,
            filter: { type: 'string' },
            renderer: symbol_list_renderer
        },{
            text: "Intron 1 Genes",
            dataIndex: 'int1sym',
            align: 'left',
            width: symWid,
            sortable: true,
            filter: { type: 'string' },
            renderer: symbol_list_renderer
        },{
            text: "Exon 2 Genes",
            dataIndex: 'ex2sym',
            align: 'left',
            width: symWid,
            sortable: true,
            filter: { type: 'string' },
            renderer: symbol_list_renderer
        },{
            text: "Intron 2 Genes",
            dataIndex: 'int2sym',
            align: 'left',
            width: symWid,
            sortable: true,
            filter: { type: 'string' },
            renderer: symbol_list_renderer
        }]
    });
    return obj;
}

function quoteattr(s, preserveCR) {
    // https://stackoverflow.com/questions/7753448/how-do-i-escape-quotes-in-html-attribute-values
    preserveCR = preserveCR ? '&#13;' : '\n';
    return ('' + s) /* Forces the conversion to string. */
        .replace(/&/g, '&amp;') /* This MUST be the 1st replacement. */
        .replace(/'/g, '&apos;') /* The 4 other predefined entities, required. */
        .replace(/"/g, '&quot;')
        .replace(/</g, '&lt;')
        .replace(/>/g, '&gt;')
    /*
      You may add other replacements here for HTML only 
      (but it's not necessary).
      Or for XML, only if the named entities are defined in its DTD.
    */ 
        .replace(/\r\n/g, preserveCR) /* Must be before the next replacement. */
        .replace(/[\r\n]/g, preserveCR);
    ;
}
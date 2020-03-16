
var start = new Date;
var bgCnt = 0;
var sosThings = {
    cntType: ["Pending","Working","Done","Failed","Unknown"],
    viewType: {
        summary: "Oligo Summary",
        overview: "Oligo Gene Overview",
        detailed: "Detailed Oligo Matches",
        setpop: "Gene Table"
    },
    setpopwin: {},
    jobdata: {},
    jobcount: 0,
    datashown: 0,
    workTime: 0,
    doneTask: {}
};
var dynamicStatus = ['date', 'type', 'meta', 'prog', 'info'];
var loadedFiles = new Object;

var cacheKeys = {
    // gene accession, symbol, description, taxa, alt symbols[]
    // human accession, human symbol, human score
    genes: [ 'g','s','d','t', 'a', 'hg', 'hs', 'sc']
};
var cachePkey = {
    genes: 'gid'
};
var objCache  = new Object;

function mmArray_renderer (val, metaData) {
    if (val == undefined) {
        metaData.tdAttr = 'bgcolor=#aaa';
        return '';
    }
    if (typeof(val) != 'object') {
        metaData.tdAttr = 'bgcolor=#a00';
        return val;//typeof(val);
    }
    return exintmm_renderer (val.length, metaData);
}

function wikipedia_renderer (val) {
    if (val == undefined) {
        return '';
    }
    return "<a href='https://en.wikipedia.org/w/index.php?search="+val+
        "' target='_blank'>"+val+"</a>";
}

function exintmm_renderer (val, metaData) {
    if (val == undefined || val == '') {
        metaData.tdAttr = 'bgcolor=#ddd';
        return null;
    }
    val = parseInt(val);
    var col = count_to_color( val );
    if (val == 0) val = null;
    metaData.tdAttr = 'bgcolor='+ col;
    return val;
}

function alias_sym_renderer (val, metaData, record) {
    var sym = record.get('s');
    if (sym == undefined) sym = record.get('sym');
    return alias_symbol_list( val, sym );
}


function count_to_color ( req ) {
    var val = parseInt( req );
    if (val == undefined|| val == 0) {
        return '#ddd';
    } else if (val == 1) {
        return'#0f0';
    } else if (val <= 5) {
        return'#cf0';
    } else if (val <= 10) {
        return'#ff0';
    } else if (val <= 50) {
        return'#f90';
    } else if (val <= 100) {
        return'#f39';
    } else if (val <= 500) {
        // #fff did not work for some reason...
        return'red';
    }
    return'red';
}

function safe_classname (txt) {
    return txt.replace(/[^a-zA-Z0-9\-_]+/g,'');
}

function exonic_renderer (val, metaData) {
    if (val == undefined) return null;
    var rv = val[0];
    if (rv == undefined) return null;
    rv = rv.replace(/\|.+/, '');
    metaData.tdCls = safe_classname( rv );
    return rv;
}

function class_as_self_renderer (val, metaData) {
    if (val == undefined) return null;
    metaData.tdCls = safe_classname( val );
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

function strand_renderer (val, metaData) {
    var cls;
    if (val == undefined) {
        cls = 'nostrand';
    } else if (val == 1) {
        cls = 'Str1';
    } else if (val == -1) {
        cls = 'Str-1';
    } else if (val == 0) {
        cls = 'Str0';
    } else {
        cls = 'StrUnk';
    }
    if (cls) metaData.tdCls = cls;
    return val;
}

function mismatch_renderer (val, metaData) {
    var cls;
    if (val == undefined) {
        cls = 'MMUnk';
    } else if (val == 0) {
        cls = 'MM0';
    } else if (val == 1) {
        cls = 'MM1';
    } else if (val == 2) {
        cls = 'MM2';
    } else {
        cls = 'MM3';
    }
    if (cls) metaData.tdCls = cls;
    return val;
}

function ambig_renderer (val, metaData) {
    var cls;
    if (val == undefined || val == 0) {
        val = '';
        cls = 'Mute';
    } else if (val >= 5) {
        cls = 'AmHi';
    } else if (val >= 3) {
        cls = 'AmMed';
    } else if (val >= 1) {
        cls = 'AmLow';
    }
    if (cls) metaData.tdCls = cls;
    return val;
}

function fancy_seq_renderer (val, metaData) {
    if (val == undefined) return null;
    var rv = "<span>";
    var bases = val.split('');
    var blen = bases.length;
    var lastCls = "";
    for (var i=0; i < blen; i++) {
        var base = bases[i];
        var clss = new Array;
        // Lower case = mismatch
        if (/[a-z]/.test(base))    clss.push('mm');
        if (/[^acgt]/i.test(base)) clss.push('am');
        var cls = clss.join(' ');
        if (cls != lastCls) {
            lastCls = cls;
            rv += "</span><span class='"+cls+"'>";
        }
        rv += base;
    }
    rv += "</span>";
    return rv;
}

function fancy_sym_renderer (val, metaData, record) {
    if (record) {
        var gid = record.get('geneid');
        if (!gid) gid = record.get('gid');
        if (gid) {
            var spl = sym_pop_link( gid );
            if (spl) return spl;
        }
    }
    return val;
}
// hg hs sc
function human_sym_renderer (val, metaData, record) {
    if (record && val) {
        var gid = record.get('geneid');
        if (!gid) gid = record.get('gid');
        if (gid) {
            var obj = gene_meta( gid );
            if (obj) return "<a gid='"+gid+"' onclick='symPop(this)' class='sym'>"+val+"</a>";

        }
    }
    return val;
}
// SAN-001791
function human_score_renderer (val, metaData, record) {
    metaData.style = 'background-color:' + perc_to_color(val);
    return val; // ? Math.floor(val) : val;
}

var scoreInflection = 80;
function perc_to_color (val) {
    var sc = parseFloat(val)
    if (sc) {
        // green - yellow - red
        var red = sc < scoreInflection ? 255 :
        Math.floor(255 * (100 - sc)/(100-scoreInflection));
        var grn = sc > scoreInflection ? 255 :
        Math.floor(255 * (sc - scoreInflection)/scoreInflection);
        return 'rgb('+red+','+grn+',0)';
    } else {
        return '#aaa';
    }
}

function locus_renderer (val, metaData, record) {
    return locus_links(val, {LL: val, GT:'GT'}) || val;
}

function locus_links (id, names, joiner) {
    var links = [];
    var defNames = { LL: 'Entrez', GT: 'GeneTracker' };
    if (!names) names = {};
    if (/^LOC\d+$/.test(id)) {
        links.push(locus_hyperlink(id, names['LL'] || defNames['LL']));
        links.push(genetracker_hyperlink(id, names['GT'] || defNames['GT']));
    }
    return links.join(joiner || ', ');
}

function locus_hyperlink (id, name) {
    if (!name) name = 'Entrez';
    return "<a href='https://www.ncbi.nlm.nih.gov/gene/?term="+id+
        "' target='_blank'>"+name+"</a>";
}

function genetracker_hyperlink (id, name) {
    if (!name) name = 'GeneTracker';
    return "<a href='http://genetracker.pri.bms.com/mev-gt/MEV-GT/#rootId/"+id+
        "' target='_blank'>"+name+"</a>";
}

function file_type_renderer (val, metaData, record) {
    if (!val) val = "";
    var jid = record.get('jid');
    if (!jid || sosThings.jobcount < 2) return val;
    var jd = get_job_data( jid );
    var num = jd.num;
    if (!num || num == undefined) num = '?';
    return val + ' #'+num;
}

function path_renderer (val, metaData, record) {
    var url = record.get('url');
    if (!url) return val || "";
    if (!val) val = "Link";
    var vHtml = val.encodeHTML();
    metaData.tdCls = 'fancy';
    var sfx = record.get('sfx') || "";
    var cls = sfx;
    var act = "target='_blank' href='"+url+"'";
    if (sfx == 'Text' || sfx == 'Fasta') {
        act += " onclick='return show_file_ajax(\""+url+"\",\""+sfx+" file\")'";
    } else if (cls == 'Params') {
        act += " onclick='return show_file_ajax(\""+url+"\",\"Param file\")'";
        var rv = "<a class='"+cls+"' "+act+">View file</a>";
        var file = record.get('path');
        if (file) rv += "&nbsp;<a class='Conf' target='_blank' href='simpleOligoSearch.pl?norun=1&jobid=0&paramfile="+file+"'>Share / Configure / Run</a>";
        return rv;
    }
    return "<a class='"+cls+"' "+act+">" + vHtml + "</a>";
}

function show_file_ajax (url, title) {
    var file = url;
    file = file.replace(/.+\//, '');
    title = title ? title + ' : ' + file : file;
    pgtm_basic_ajax_launch(url, { isajax: 1 }, function (rsp) {
        var txt = text_from_response(rsp);
        raw_text_window( txt, title);
    });
    return false;
}

function text_from_response (rsp) {
    if (!rsp) return "<i class='CodeError'>No response from server!</i>";
    var txt = rsp.responseText;
    if (txt) return txt;
    var stat = rsp.status;
    if (stat == 200) {
        return "<i class='note'>File is empty</i>";
    }
    txt = "<i class='warn'>No data recovered from file\n";
    txt += "Status: " +rsp.statusText + "\n";
    txt += "Code: " + rsp.status;
    txt += "</i>";
    return txt;
}

var mmWid = 40;
var symWid = 200;
var colKeyAliases = {
    s:   'sym',
    t:   'taxa',
    gid: 'geneid',
    g:   'geneacc',
    d:   'genedesc'
};
var coltypes = {
    oligoid: {
        text: "OligoID",
        dataIndex: 'oligoid',
        width: 100,
        sortable: true,
        cellWrap: false,
        filter : { type : 'list' }
    },
    id: {
        text: "DB ID",
        dataIndex: 'id',
        width: 100,
        sortable: true,
        cellWrap: false,
        hidden: true
    },
    oligoseq : {
        text: "OligoSeq",
        dataIndex: 'oligoseq',
        align: 'left',
        width: 150,
        tdCls: 'Seq',
        sortable: true,
        filter: { type: 'string' }
    },
    subseq : {
        text: "GenomicSeq",
        dataIndex: 'subseq',
        align: 'left',
        width: 150,
        tdCls: 'Seq',
        renderer: fancy_seq_renderer,
        sortable: true,
        filter: { type: 'string' }
    },
    sym: {
        text: "Symbol",
        dataIndex: 'sym',
        width: 100,
        renderer: fancy_sym_renderer,
        sortable: true,
        filter : { type : 'string' }
    },
    a: {
        text: "Symbol Aliases",
        dataIndex: 'a',
        width: 300,
        sortable: false,
        hidden: true,
        cls: 'altsym',
        renderer: alias_sym_renderer,
        filter : { type : 'string' }
    },
    mm: {
        text: "MM",
        dataIndex: 'mm',
        align: 'center',
        width: 40,
        sortable: true,
        filter: { type: 'list' },
        renderer: mismatch_renderer
    },
    str: {
        text: "Str",
        dataIndex: 'str',
        align: 'center',
        width: 40,
        sortable: true,
        filter: { type: 'list' },
        renderer: strand_renderer
    },
    exon: {
        text: "Exonic",
        dataIndex: 'exon',
        align: 'center',
        width: 100,
        sortable: true,
        filter: { type: 'list' },
        renderer: exonic_renderer
    },
    notes: {
        text: "!!",
        dataIndex: 'notes',
        align: 'center',
        width: 30,
        sortable: true,
        filter: { type: 'string' },
        getTip: function () { return "hello"; },
        renderer: note_renderer
    },
    alignment: {
        text: "Alignment",
        dataIndex: 'alignment',
        align: 'left',
        width: 200,
        tdCls: 'Seq',
        sortable: false
    },
    genedesc: {
        text: "Gene Description",
        dataIndex: 'genedesc',
        align: 'left',
        width: 300,
        sortable: true,
        filter: { type: 'string' }
    },
    taxa: {
        text: "Species",
        dataIndex: 'taxa',
        align: 'left',
        width: 200,
        sortable: true,
        tdCls: 'tax',
        renderer: wikipedia_renderer,
        hidden: false
    },
    geneid: {
        text: "GeneId",
        dataIndex: 'geneid',
        align: 'center',
        width: 100,
        sortable: true,
        hidden: true
    },
    geneacc: {
        text: "Accession",
        dataIndex: 'geneacc',
        align: 'left',
        width: 100,
        renderer: locus_renderer,
        sortable: true,
        hidden: true
    },
    humanacc: {
        text: "Human Acc",
        dataIndex: 'humanacc',
        align: 'left',
        width: 100,
        renderer: locus_renderer,
        sortable: true
    },
    humansym: {
        text: "Hs Sym",
        dataIndex: 'humansym',
        width: 100,
        renderer: human_sym_renderer,
        sortable: true,
        filter : { type : 'string' }
    },
    humanscore: {
        text: "Hs %",
        dataIndex: 'humanscore',
        width: 40,
        renderer: human_score_renderer,
        sortable: true
    },
    footprint: {
        text: "Genomic Footprint",
        dataIndex: 'footprint',
        align: 'left',
        width: 200,
        sortable: true,
        hidden: true,
        filter: { type: 'string' }
    },
    am: {
        text: "AM",
        dataIndex: 'am',
        align: 'center',
        width: 40,
        sortable: true,
        hidden: false,
        filter: { type: 'list' },
        renderer: ambig_renderer
    },
    rnafoot: {
        text: "RNA Footprint",
        dataIndex: 'rnafoot',
        align: 'left',
        width: 200,
        sortable: true,
        hidden: true,
        filter: { type: 'string' },
        renderer: rnafoot_renderer
    },
    exonnum: {
        text: "Exon#",
        dataIndex: 'exonnum',
        align: 'center',
        width: 90,
        sortable: true,
        hidden: true,
        filter: { type: 'string' }
    },
    rnanum: {
        text: "#RNAs",
        dataIndex: 'rnanum',
        align: 'center',
        width: 60,
        sortable: true,
        hidden: false,
        filter: { type: 'numeric' }
    },
    sbjid: {
        text: "Subject ID",
        dataIndex: 'sbjid',
        align: 'left',
        width: 200,
        sortable: true,
        hidden: true
    },
    sbjpos: {
        text: "Subject Pos",
        dataIndex: 'sbjpos',
        align: 'left',
        width: 80,
        sortable: true,
        hidden: true
    },
    build: {
        text: "Build",
        dataIndex: 'build',
        align: 'center',
        width: 100,
        sortable: true,
        renderer: class_as_self_renderer,
        filter: { type: 'list' }
    },
    ex0sym: {
        text: "Exon Perfect",
        dataIndex: 'ex0sym',
        align: 'left',
        torcls: 'MM0',
        width: symWid,
        sortable: false,
        renderer: symbol_list_renderer
    },
    ex1sym: {
        text: "Exon 1MM",
        dataIndex: 'ex1sym',
        align: 'left',
        torcls: 'MM1',
        width: symWid,
        sortable: false,
        renderer: symbol_list_renderer
    },
    ex2sym: {
        text: "Exon 2MM",
        dataIndex: 'ex2sym',
        align: 'left',
        torcls: 'MM2',
        width: symWid,
        sortable: false,
        renderer: symbol_list_renderer
    },
    int0sym: {
        text: "Intron Perfect",
        dataIndex: 'int0sym',
        torcls: 'MM0',
        align: 'left',
        width: symWid,
        sortable: false,
        renderer: symbol_list_renderer
    },
    int1sym: {
        text: "Intron 1MM",
        dataIndex: 'int1sym',
        align: 'left',
        torcls: 'MM1',
        width: symWid,
        sortable: false,
        renderer: symbol_list_renderer
    },
    int2sym: {
        text: "Intron 2MM",
        dataIndex: 'int2sym',
        align: 'left',
        torcls: 'MM2',
        width: symWid,
        sortable: false,
        renderer: symbol_list_renderer
    },
    ex0: {
        text: "Ex0",
        dataIndex: 'ex0',
        align: 'center',
        width: mmWid,
        sortable: true,
        filter: { type: 'numeric' },
        renderer: exintmm_renderer
    },
    ex1: {
        text: "Ex1",
        dataIndex: 'ex1',
        align: 'center',
        width: mmWid,
        sortable: true,
        filter: { type: 'numeric' },
        renderer: exintmm_renderer
    },
    ex2: {
        text: "Ex2",
        dataIndex: 'ex2',
        align: 'center',
        width: mmWid,
        sortable: true,
        filter: { type: 'numeric' },
        renderer: exintmm_renderer
    },
    int0: {
        text: "Int0",
        dataIndex: 'int0',
        align: 'center',
        width: mmWid,
        sortable: true,
        filter: { type: 'numeric' },
        renderer: exintmm_renderer
    },
    int1: {
        text: "Int1",
        dataIndex: 'int1',
        align: 'center',
        width: mmWid,
        sortable: true,
        filter: { type: 'numeric' },
        renderer: exintmm_renderer
    },
    int2: {
        text: "Int2",
        dataIndex: 'int2',
        align: 'center',
        width: mmWid,
        sortable: true,
        filter: { type: 'numeric' },
        renderer: exintmm_renderer
    },
    mm0: {
        text: "MM0",
        dataIndex: 'mm0',
        align: 'center',
        width: mmWid,
        sortable: true,
        filter: { type: 'numeric' },
        renderer: exintmm_renderer
    },
    mm1: {
        text: "MM1",
        dataIndex: 'mm1',
        align: 'center',
        width: mmWid + 5,
        sortable: true,
        filter: { type: 'numeric' },
        renderer: exintmm_renderer
    },
    mm2: {
        text: "MM2",
        dataIndex: 'mm2',
        align: 'center',
        width: mmWid + 10,
        sortable: true,
        filter: { type: 'numeric' },
        renderer: exintmm_renderer
    },
    database: {
        text: "Database Path",
        dataIndex: 'database',
        align: 'left',
        width: 300,
        sortable: true,
        hidden: true,
        filter: { type: 'string' }
    }
};

for (var cka in colKeyAliases) {
    var src = coltypes[ colKeyAliases[cka] ];
    var data = coltypes[cka] = new Object;
    for (var k in src) {
        data[k] = src[k];
    }
    data.dataIndex = cka;
}

/* Lots of ExtJS stuff is being guided by code written by Mark Russo,
   in some cases taken in whole cloth from it

   http://kraken.pri.bms.com/biohtml/russom/dev/aso/public_html/extjsgrid2.html

   ExtJS / Sencha docs:
   Store:    http://docs.sencha.com/extjs/4.2.1/#!/api/Ext.data.Store
   MemProxy: http://docs.sencha.com/extjs/4.2.1/#!/api/Ext.data.proxy.Memory
   
*/

// Ext.QuickTips.init();

// https://stackoverflow.com/questions/7918868/how-to-escape-xml-entities-in-javascript
if (!String.prototype.encodeHTML) {
    String.prototype.encodeHTML = function () {
        return this.replace(/&/g, '&amp;')
            .replace(/</g, '&lt;')
            .replace(/>/g, '&gt;')
            .replace(/"/g, '&quot;')  //" 
            .replace(/'/g, '&apos;'); //'
    };
}

function sos_error (txt) {
    if (!txt) return;
    alert(txt);
    return;
    var panel = error_panel();
    var msg   = newExt.dom.Element( 'pre' );
    msg.setHTML( txt );
    errPanel.addChildEls( msg );
}

function monitor_pgtm (jid) {
    // Initiates background monitoring of a particular job
    var pgtm = new PgTaskManager
    ( process_ss_json, { "jobid": jid, "getstatus": 1 } );
    pgtm.error_callback( sos_error );
}

function process_ss_json ( data, pgtm, url, params ) {
    //if (!pgtm) {
    //    if (!data) data = {};
    //    data.jsonError = "No PGTM object defined";
    //    dump_json_to_status( data );
    //    return null;
    //}
    if (!data) {
        var err = "Empty server response!"
        data = { jsonError: err };
        if (pgtm) {
            pgtm.error_msg(err + " : " + pgtm.paramsToGetString());
        } else {
            dump_json_to_status( data );
        }
        return;
    }
    var type = data.type;
    cache_data( data );
    if (!type) {
        var err = "Server response does not specify data type";
        if (pgtm) {
            pgtm.error_msg(err + " : " + pgtm.paramsToGetString());
        } else {
            data.jsonError = err;
            dump_json_to_status( data );
        }
        return;
    }
    if (type == 'jobstatus') {
        return process_pgtm_status( data, pgtm );
    } else if (sosThings.viewType[ type ]) {
        return add_to_view( data, pgtm, url, params );
    }
    if (pgtm) pgtm.error_msg("Unrecognized server response type '"+type+"' : "+
                             pgtm.paramsToGetString());
    return null;
}

function cache_data ( data ) {
    if (!data) return;
    for (var key in cacheKeys) {
        var dc  = data[ key ];
        if (!dc) continue;
        // The data contains supporting information
        var tc = objCache[ key ];
        if (!tc) tc = objCache[ key ] = new Object;
        var subKeys = cacheKeys[key];
        var skl     = subKeys.length;
        var pk      = cachePkey[key];
        // Copy each supporting data over to local cache
        for (var id in dc) {
            if (tc[id]) continue; // Already cached
            // Add a new object to the cache
            var trg = tc[ id ] = { id: id };
            trg[ pk ] = id;
            var src = dc[ id ];
            for (var i = 0; i < skl; i++) {
                var v = src[i];
                trg[ subKeys[i] ] = (v == undefined) ? "" : v;
            }
        }
    }
}

function gene_meta (gid) {
    if (!gid) return 0;
    var tc = objCache.genes || new Object;
    return (gid in tc) ? tc[gid] : null;
}

function gene_symbol (gid) {
    var obj = gene_meta( gid );
    if (obj) {
        return obj.s || obj.g || 'NoACC!';
    } else if (obj == undefined) {
        return 'Unk!';
    }
    return "";
}

function gene_desc ( gid ) {
    var obj = gene_meta( gid );
    if (obj) {
        return obj.d || "";
    } else if (obj == undefined) {
        return 'Unkown gene_id '+gid+' !';
    }
    return "";
}

function note_file( data ) {
    if (!data) return;
    fs = file_store();
    fp = file_panel();
    var path = data.path;
    if (path == undefined) path = data.url;
    if (path) {
        if (/\.txt$/.test(path)) {
            data.sfx = 'Text';
        } else if (/\.fa$/.test(path)) {
            data.sfx = 'Fasta';
        } else if (/\.tsv$/.test(path)) {
            data.sfx = 'TSV';
        } else if (/\/$/.test(path)) {
            data.sfx = 'Folder';
        } else if (/\.xlsx?$/.test(path)) {
            data.sfx = 'Excel';
        } else if (/\.param$/.test(path)) {
            data.sfx = 'Params';
        }
    }
    // if (!data.sfx) dumper_window( data );
    fs.add( data );
    fp.getView().refresh();
    // fp.doLayout();
}

function status_renderer (val, metaData, record) {
    // .type = simple Done/Working/Pending etc
    // .prog = status column from pgtm_status
    if (val == undefined) return "?";
    var prog = record.get('prog');
    if (!prog) {
        // No progress data
        metaData.tdCls = 'stat' + val;
        return val;
    }
    var subtask = prog.subtask;
    if (!subtask) return val;
    var i = subtask.length - 1;
    while (i > 0) {
        if (subtask[i][1] != 0) break;
        i--;
    }
    var st  = subtask[i];
    var txt = st[0];
    if (st[2]) txt += ':'+st[2];
    //metaData.tdCls = 'progbar';
    var perc = Math.floor( (st[1] || 0) * 100);
    var rv = "<div class='progbar'><span class='statWorking'>"+txt+"</span><div style='width:"+perc+"%'></div></div>";
    return rv;
}

function process_pgtm_status ( data, pgtm ) {
    var statDiv = status_div();
    var ss = status_store();
    var sp = status_panel();
    var title = new Object;
    var tott  = data.count.Total;
    var jid   = data.jid;
    if (data.err) {
        ss.add( { type: "Failed",
                  tok: "",
                  status: data.err } );
    }
    sosThings.cntType.map( function (v) {
        title[ v ] = 0;
    });
    
    // Note tasks already stored
    var oldRecs = new Object;
    ss.each( function (rec) {
        var tid = rec.get('tid');
        if (!tid) return;
        oldRecs[ tid ] = rec;
        var ojid = rec.get('jid');
        if (ojid && ojid != jid) {
            // Note progress from other jobs
            var otype = rec.get('type') || 'Unknown';
            title[ otype ]++;
        }
    });
    // dumper_window( data, 'PGTM Status' );
    var stats = data.status || [];
    var snum  = stats.length;
    for (var i=0; i < snum; i++) {
        var rec  = stats[i];
        var tid  = rec.tid;
        var old  = oldRecs[tid];
        var type = rec.type || "Unknown";
        title[type]++;
        if (old) {
            // The task is already in the store
            // If nothing has changed, nothing to do:
            if (old.get('date') == rec.date) continue;
            // Update the old entry
            for (var j = 0; j < dynamicStatus.length; j++) {
                var field = dynamicStatus[ j ];
                old.set(field, rec[ field ]);
            }
        } else {
            // New task, add to store
            ss.add( rec );
        }
        if (type == 'Done') {
            // We may have output that needs to be handled
            var url = rec.url;
            if (url) {
                if (! sosThings.doneTask[tid]++) {
                    pgtm.launch_ajax( url, {} );
                }
            }
        }
    }
    // Make the title
    var tiBits = new Array();
    var ct     = sosThings.cntType.length;
    for (var i = 0; i < ct; i++) {
        var type = sosThings.cntType[ i ];
        var num  = title[type];
        if (!num) continue;
        var cls  = 'stat'+type;
        tiBits.push("<span class='"+cls+"'>"+num+"</span>");
    }
    var tiTxt = tiBits.join('/');
    var todo = title.Working + title.Pending;
    if (todo || data.wait) {
        // Some tasks have not finished, or all tasks are done but we
        // need to wait for the parent process to assemble them
        var why = todo ? todo + ' jobs running' : data.wait;
        query_alert('Job '+jid+': '+why);
        // check back in a second:
        setTimeout(function () {
            pgtm.launch_ajax( );
        }, 1000 );
    } else {
        // All appear to have completed
        tiTxt += " <i>Done</i>";
        // If the status tab is active switch to the content tab
        query_alert('Job '+jid+': Complete');
    }
    var now = new Date;
    var min = Math.floor(0.5 + (now.getTime() - start.getTime()) / 600) / 100;
    tiTxt += " - " + min + "min";
    sp.setTitle(tiTxt);
    // sp.doLayout();

    if (data.path) {
        note_file( { type: "Results", jid: jid,
                     path: data.path, url: data.url } );
        if (data.url) window.open( data.url, '_blank' );
    }
    sp.getView().refresh();
    // dump_json_to_status( data );
}

function dump_json_to_status ( data ) {
    var statDiv = status_div();
    var txt = JSON.stringify( data );
    statDiv.innerHTML = "<pre>" + txt.encodeHTML() + "</pre>";
}

function status_div () {
    if (!sosThings.statusDiv) {
        sosThings.statusDiv = document.createElement('div');
        document.body.appendChild( sosThings.statusDiv );
    }
    return sosThings.statusDiv;
}

function status_store () {
    if (!sosThings.statusStore) {
        sosThings.statusStore = Ext.create('Ext.data.Store', {
            storeId:'statusStore',
            fields:['type', 'tok', 'status'],
            data:{'items':[ ]},
            proxy: {
                type: 'memory',
                reader: {
                    type: 'json',
                    rootProperty: 'items'
                }
            }
        });
    }
    return sosThings.statusStore;
}

function status_panel () {
    var panel = sosThings.statusPanel;
    if (panel) return panel;
    var tf      = tab_frame();
    var store   = status_store();
    panel       = sosThings.statusPanel = Ext.create('Ext.grid.Panel', {
        title: 'Analysis Status',
        width: 910,
        height: 100,
        store: store,
        columns: [
            { text: 'Status',
              width: 150,
              renderer: status_renderer,
              dataIndex: 'type'
            },
            { text: 'Task', width: 250, dataIndex: 'tok' },
            { text: 'Information', width: 500, dataIndex: 'info' },
            { text: 'JobID', width: 100, dataIndex: 'jid', hidden: true },
            { text: 'TaskD', width: 100, dataIndex: 'jid', hidden: true }
       ]
    });
    tf.add( panel );
    // tf.setActiveTab( panel );
    return panel;
}

function query_alert (msg) {
    var el = document.getElementById('srchmsg');
    if (el) {
        el.innerHTML = msg;
    } else {
        alert(msg);
    }
}

function handle_query_response (rsp, url, params) {
    var data = pgtm_basic_json_parse (rsp, url, params, sos_error);
    if (!data) return;
    if (data.alert) query_alert(data.alert);
    if (data.jobid) {
        monitor_pgtm( data.jobid );
        set_job_data( data );
    }
}

function set_job_data (data) {
    if (!data) return;
    var jid = data.jobid;
    if (!jid) return;
    var jd = sosThings.jobdata[ jid ] = data;
    jd.num = ++sosThings.jobcount;
    if (data.files) {
        var urls = data.urls || {};
        for (var key in data.files) {
            var url = urls[ key ];
            var path = data.files[ key ];
            note_file( { type: key,
                         jid: jid,
                         url: url,
                         path: path } );
        }
    }
}

function get_job_data (jid) {
    return sosThings.jobdata[ jid ] || {};
}

function query_panel (domid) {
    var panel = sosThings.queryPanel;
    if (panel) return panel;
    var form = document.getElementById(domid);
    if (!form) return;
    var tf    = tab_frame();
    var panel = sosThings.queryPanel = Ext.create('Ext.panel.Panel', {
        title: 'Query',
        contentEl: domid
    });
    tf.add( panel );
    var url = window.location.pathname;
    form.onsubmit = function () {
        var params = mini.form.gethash(form);
        params.dhtmlajax = 1;
        pgtm_basic_ajax_launch(url, params, handle_query_response);
        query_alert('Request submitted ...');
        return false;
    };
    return panel;
}

function error_panel () {
    var panel = sosThings.errorPanel;
    if (panel) return panel;
    if (!domid) return;
    var tf    = tab_frame();
    var panel = sosThings.errorPanel = tf.add({
        xtype: 'box',
        cls: 'ajaxTaskErrors',
        autoEl: {cn: text}
    });
    panel.setTitle('Errors');
    return panel;
}

function file_store () {
    if (!sosThings.fileStore) {
        sosThings.fileStore = Ext.create('Ext.data.Store', {
            storeId:'fileStore',
            fields:[ 'type', 'jid', 'path', 'url' ],
            data:{'items':[ ]},
            proxy: {
                type: 'memory',
                reader: {
                    type: 'json',
                    rootProperty: 'items'
                }
            }
        });
    }
    return sosThings.fileStore;
}

function file_panel () {
    var panel = sosThings.filePanel;
    if (panel) return panel;
    var tf      = tab_frame();
    var store   = file_store();
    panel       = sosThings.filePanel = Ext.create('Ext.grid.Panel', {
        title: 'Files',
        width: 910,
        height: 100,
        store: store,
        plugins: 'gridfilters', // adds a simple and flexible filter set
        columns: [
            { text: 'JobID',
              width: 50,
              filter : { type : 'list' },
              dataIndex: 'jid' },
            { text: 'Suffix',
              width: 60,
              filter : { type : 'list' },
              dataIndex: 'sfx' },
            { text: 'Type',
              width: 150,
              renderer: file_type_renderer,
              filter : { type : 'list' },
              dataIndex: 'type' },
            { text: 'Location',
              width: 400,
              renderer: path_renderer,
              dataIndex: 'path'
            },
            { text: 'Hyperlink', hidden: true, width: 300, dataIndex: 'url' }
        ]
    });
    tf.add( panel );
    // tf.setActiveTab( panel );
    return panel;
}

function add_to_view (data, pgtm, url, params) {
    if (!data) return;
    if (!sosThings.doneTask[ url ]) sosThings.doneTask[ url ] = 0;
    if (sosThings.doneTask[ url ]++) return;
    // pgtm.error_msg(sosThings.doneTask[ url ] + " = '"+ url+"'" );
    var sp = store_panel(data);
}

function store_panel( data ) {
    if (!data) return [];
    var type = data.type || "Custom Data";
    var store = store_for_type( type );
    store.add( data.rows );
    var panel = panel_for_type( type, data );
    // panel.doLayout();
    set_panel_title(panel)
    panel.getView().refresh();
    return [store, panel];
}

function set_panel_title(panel) {
    if (!panel) return;
    var store = panel.store;
    var type  = panel.CATtype
    var vname = sosThings.viewType[ type ] || type;
    var snum  = store.data.getCount();
    var vnum  = store.getCount();
    vname += ' - ' + snum + " row" + (snum == 1 ? '' : 's');
    if (vnum < snum) {
        vname += ', '+vnum+' match filter';
    }
    panel.setTitle( vname );    
}

function store_for_type (type, sorters) {
    var sid   = "STORE:"+type;
    var store = sosThings[sid];
    if (store) return store;
    store = sosThings[sid] = Ext.create('Ext.data.Store', {
        storeId: sid,
        fields:[],
        data:{'items':[ ]},
        buffered: false,
        sorters: sorters,
        proxy: {
            type: 'memory',
            reader: {
                type: 'json',
                rootProperty: 'items'
            }
        }
    });
    return store;
}

function normalized_columns ( srcCols, oldCols ) {
    // srcCols should be an array of strings corresponding to ordered
    // data indices
    // oldCols is an optional array of prior column hashes
    var hasCol  = new Object;
    if (oldCols) {
        for (var j = 0; j < oldCols.length; j++) {
            if (oldCols[j] && oldCols[j].dataIndex) {
                hasCol[ oldCols[j].dataIndex ] = 1;
            }
        }
    }
    var newCols = new Array;
    for (var c=0; c < srcCols.length; c++) {
        var cn = srcCols[c];
        if (hasCol[ cn ]) {
            // This column is already present
            continue;
        }
        var cInfo = coltypes[cn];
        if (!cInfo) {
            var name = cn;
            name = name.replace('USER-','');
            cInfo = {
                text: name,
                width: 150,
                dataIndex: cn
            };
        }
        newCols.push(cInfo);
    }
    return newCols;
}

function tab_frame () {
    var tf = sosThings.tabFrame;
    if (tf) return tf;
    var targ = 'extout';
    tf = sosThings.tabFrame = Ext.create('Ext.tab.Panel', {
        renderTo: targ,
        width: 1000,
        height: 600
    });
    return tf;
}

function panel_for_type ( type, data ) {
    var pid   = "PANEL:"+type;
    var oldPanel = sosThings[pid];
    if (!data) return oldPanel;
    var oldCols = oldPanel ? oldPanel.columns : null;

    var newCols = normalized_columns( data.cols, oldCols );
    
    if (newCols.length == 0) {
        // No new columns
        if (!oldPanel) alert("No columns defined for "+type+" data");
        return oldPanel;
    }
    var cols = newCols;
    if (oldCols) {
        // A panel already existed, but we now have new columns
        // we need to make a new panel with extended columns
        cols = oldCols.concat(newCols);
        oldPanel.ownerCt.remove( oldPanel );
    }
    var tf    = tab_frame();
    var store = store_for_type( type );
    var panel = sosThings[pid] = new Ext.grid.Panel({
        title: 'Results',
        width: 1000,             // Dimensions
        height: 500,
        // renderTo: tf,
        collapsible: true,      // collapes using double arrows in header.
        store: store,           // The store object that manages data.
        multiSelect: true,      // Allows the user to select multiple rows
        plugins: 'gridfilters', // adds a simple and flexible filter set
        
        viewConfig: {
            loadMask: true, //false
            loadingText: 'Downloading results...'
        },
        selModel: {
            pruneRemoved: false
        },
        columns: cols
    });
    tf.add(panel);
    panel.CATtype = type;
    return panel;
}

var maxSymsToShow = 3;

function symbol_list_renderer  (val, metaData, record, ri, ci) {
    if (val == undefined || !val) return '';
    var to = typeof(val);
    if (to != 'object') return '?';
    var list = new Array;
    var lnum = val.length;
    var max  = lnum > maxSymsToShow ? maxSymsToShow : lnum;
    for (var i = 0; i < max; i++) {
        var spl = sym_pop_link( val[i] );
        if (spl) list.push( spl );
    }
    var rv = list.join(', ');
    var col = count_to_color( lnum );
    if (lnum > 1) {
        // Multiple entries here
        var rid = record.getId();
        rv = "<a rid='"+rid+"' ci='"+ci+"' onclick=setPop(this) style='background-color:"+col+"' class='setpop' title='Click to see all genes'>"+lnum+"</a> "+rv;
        if (lnum > maxSymsToShow) rv += '&hellip;'
    }
    // metaData.tdAttr = 'bgcolor='+ col;
    return rv;
}

function sym_pop_link (gid) {
    if (!gid) return "";
    var obj = gene_meta( gid );
    if (!obj) return "<span class='CodeError'>No gene_id="+gid+"</span>";
    var sym = obj.s;
    var cls = 'sym';
    if (!sym) {
        // Use the accession if no symbol is available
        sym = obj.g;
        cls = 'symless';
        if (!sym) {
            // Oopsies
            sym = 'NoACC!';
            cls = 'CodeError';
        }
    }
    var title = gene_desc( gid );
    if (title) title = " title='"+quoteattr(title)+"'";
    return "<a gid='"+gid+"' onclick='symPop(this)' class='"+cls+"'"+
        title+">"+sym+"</a>";
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

var spCols = ['gid', 's', 'g', 't', 'd', 'a' ];
function setPop (anc) {
    if (!anc) return false;
    var rid = anc.getAttribute('rid');
    var ci  = anc.getAttribute('ci');
    if (rid == undefined || ci == undefined) return false;
    var win = sosThings.setpopwin[ rid ];
    if (win) {
        win.toFront();
        return;
    }
    
    var os = store_for_type( 'overview' );
    if (!os) return false;
    var record = os.getById( rid );
    if (!record) return false;
    var op = panel_for_type( 'overview' );
    if (!op) return false;
    var colD = op.columns[ ci ];
    var cn   = colD.dataIndex;
    var list = record.get(cn);
    if (!list) {
        alert("Failed to recover gene list for row:"+rid+" col:"+cn);
        return false;
    }
    var llen = list.length;
    if (!llen) {
        alert("Gene list is empty for row:"+rid+" col:"+cn);
        return false;
    }
    var data = new Array;
    for (var i = 0; i < llen; i++) {
        var gm = gene_meta( list[i] );
        if (gm) data.push(gm);
    }
    var type  = 'setpop';
    var store = store_for_type( type + rid, [{property: 's'}] );
    store.removeAll();
    store.add( data );
    var cols  = normalized_columns( spCols );
    var tbits = new Array;
    if (sosThings.viewType[ type ]) tbits.push(sosThings.viewType[ type ]);
    var oid = record.get('oligoid');
    if (oid) tbits.push(oid);
    if (colD.text) {
        var tcls = safe_classname( colD.torcls ) || "";
        tbits.push( "<span class='tor "+tcls+"'>"+colD.text+"</span>" );
    }
    var bld = record.get('build');
    if (bld) {
        var tcls = safe_classname( bld );
        tbits.push( "<span class='tor "+tcls+"'>"+bld+"</span>" );
    }
    tbits.push( llen + ' entr'+(llen == 1 ? 'y' : 'ies') );
    var title = tbits.join(' - ');
    
    win = sosThings.setpopwin[ rid ] = Ext.create('Ext.window.Window', {
        title: title,
        height: 500,
        width: 600,
        layout: 'fit',
        items:[{
            xtype : 'grid',
            autoScroll: true,
            plugins: 'gridfilters',
            columns: cols,
            store: store

        }]
    }).show();
}

function symPop (anc) {
    if (!anc) return false;
    var gid = anc.getAttribute('gid');
    if (!gid) return false;
    var title = anc.innerHTML;
    if (anc.title) title = title + " : " + anc.title;
    var html = "";
    var obj = gene_meta( gid );
    if (!obj) {
        html = "?? Could not locate data for gene_id = "+gid+" ??";
    } else {
        var acc  = obj.g;
        if (acc) {
            html += "<h5>"+acc+"</h5>";
            html += "<b>Links:</b> "+locus_links(acc)+"<br />\n";
        } else {
            html += "<span class='alert'>Accession not found!</span><br />\n";
        }
        var sym = obj.s;
        if (sym) {
            html += "<b>Symbol:</b> <span class='sym'>"+sym+"</span>";
            var ali = alias_symbol_list( obj.a, sym );
            if (ali) html += "<span class='altsym'>, " +ali+ "</span>";
            html += "<br />\n";
        }
        var desc = obj.d;
        if (desc) {
            html += "<b>Description:</b> <span class='desc'>"+desc+"</span><br />\n";
        }
        var tax = obj.t;
        if (tax) {
            html += "<b>Taxa:</b> <span class='tax'>"+tax+"</span><br />\n";
        }
        var hg = obj.hg
        if (hg && hg != obj.g) {
            var perc = obj.sc;
            html += "<b>Human Orthologue:</b> "+hg+" <span style='background-color:"+perc_to_color(perc)+"'>"+perc+"%</span>";
            if (obj.hs) html += " <i>"+obj.hs+"</i>";
            html += ' ' + locus_links(hg)
            html += "<br />\n";
        }
    }
    // html += "<pre>" + JSON.stringify( obj, null, 2 ) + "</pre>";
    var win = Ext.create('Ext.window.Window', {
        title: title,
        height: 500,
        width: 600,    
        items:[{
            xtype: 'box',
            html: html
        }]
    }).show();
    return true;
}

function dumper_window (obj, subtitle) {
    // if (!obj) return;
    var win = sosThings.dumpWindow;
    if (!win) {
        win = sosThings.dumpWindow = Ext.create('Ext.window.Window', {
            autoScroll: true,
            listeners: {
                beforeclose : function () {
                    sosThings.dumpWindow = null;
                }
            },
            height: 500,
            width: 600
        });
        win.otherHits = 0;
    }
    if (win.isVisible()) {
        win.otherHits++;
        return;
    }
    var title = 'Object Debuger';
    if (subtitle) title += ' : ' + subtitle;
    win.setTitle(title);
    var html = "<pre>" + JSON.stringify( obj, null, 2 ) + "</pre>";
    win.removeAll();
    win.add( [{ xtype: 'box', html: html }] );
    win.show();
}

function raw_text_window (txt, title) {
    var win = Ext.create('Ext.window.Window', {
        autoScroll: true,
        height: 500,
        width: 600,
        items:[{
            xtype: 'box',
            html: "<pre>"+txt+"<pre>"
        }]
    });
    if (!title) title = "Text Dump Window";
    win.setTitle(title);
    win.show();
}

function alias_symbol_list (ali, sym) {
    if (!ali) return "";
    var list = new Array;
    for (var i = 0; i < ali.length; i++) {
        var a = ali[i];
        if (a && a != sym) list.push(a);
    }
    return list.sort().join(', ');
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

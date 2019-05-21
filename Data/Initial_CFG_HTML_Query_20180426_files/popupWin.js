var popupWin

Stamp = new Date();
var Hours;
var Mins;
var Seconds;
var uniqueWindow;
Hours = Stamp.getHours();
Minutes = Stamp.getMinutes();
Seconds = Stamp.getSeconds();
uniqueWindow = "popupWin_" + Hours + Minutes + Seconds;
var win;

//Function for opening a new window of given width or height
function openWindow(page, wid, hei) {
    if (wid == null) {
        wid = 1200;
    }
    if (hei == null) {
        hei = 600;
    }
    if (!popupWin || popupWin.closed) {
        popupWin = window.open(page, uniqueWindow, "width=" + wid + ",height=" + hei + ",resizable=yes,scrollbars=yes");
        popupWin.moveTo(250, 250);
        popupWin.focus();
    } else {
        /*popupWin.close(); */
        popupWin = window.open(page, uniqueWindow, "width=" + wid + ",height=" + hei + ",resizable=yes,scrollbars=yes");
        popupWin.focus();
    }
}

//Opening a new window
function showPopup(some_url)
{
    win = window.open(some_url, 'popup', 'width=500,height=220,scrollbars=no,menubar=no');
    win.moveTo(x + 20, y - 5);
    return true;
}

//Closing a popup window.
function closePopup()
{
    win.close();
    return true;
}
//end hide --></script></head>

//
function closeMe(url) {
    var wopener = url;

	if(window.opener)  {
        window.opener.location.href = wopener;
    	window.close();
    } else {
    	window.location.href = wopener;
    }
}
function move_box(an, box) {
    var cleft = 0;
    var ctop = 0;
    var obj = an;
    while (obj.offsetParent) {
        cleft += obj.offsetLeft;
        ctop += obj.offsetTop;
        obj = obj.offsetParent;
    }
    //box.style.left = cleft + 'px';
    box.style.left = 0 + 'px';
    ctop += an.offsetHeight + 8;
    if (document.body.currentStyle &&
        document.body.currentStyle['marginTop']) {
        ctop += parseInt(
                document.body.currentStyle['marginTop']);
    }
    //box.style.top = ctop + 'px';
    box.style.top = 0 + 'px';
}
//
//function move_box( box, left, top) {
//
//    box.style.left = x + left + 'px';
//	box.style.top = y + top + 'px';
//}

//
function fixed_box(box, left, top)
{
    //box.style.left = left + 'px';
    //box.style.top = top + 'px';
    box.style.left = 0 + 'px';
    box.style.top = 0 + 'px';
}

function show_hide_box(an, width, height, borderStyle) {
    hide("br");
    var href = an.href;
    var boxdiv = document.getElementById("br");

    if (boxdiv != null) {
        if (boxdiv.style.display == 'none') {
            move_box(an, boxdiv, 0, 0);
            boxdiv.style.display = 'block';
        } else
            boxdiv.style.display = 'none';
        return false;
    }

    boxdiv = document.createElement('div');
    boxdiv.setAttribute('id', "br");
    // "br"

    boxdiv.style.display = 'block';
    boxdiv.style.position = 'absolute';
    boxdiv.style.width = width + 'px';
    boxdiv.style.height = height + 'px';
    boxdiv.style.border = borderStyle;
    boxdiv.style.backgroundColor = '#fff';

    var contents = document.createElement('iframe');


    contents.scrolling = 'no';
    contents.frameBorder = '0';
    contents.style.width = width + 'px';
    contents.style.height = height + 'px';
    contents.style.backgroundColor = "#ccc"
    contents.src = href;

    boxdiv.appendChild(contents);

    document.body.appendChild(boxdiv);
    move_box(an, boxdiv, 0, 0);

    return false;
}

//
function hide(an) {
    //var href = an;
    var boxdiv = document.getElementById(an);
    if (boxdiv != null) {
        if (boxdiv.style.display == 'block')
            boxdiv.style.display = 'none';
    }
}

function hide2() {
    //var href = an;
    var boxdiv = document.getElementById("br");
    var select = document.getElementById("select")

    if (boxdiv.style.display == 'block')
        boxdiv.style.display = 'none';

    if (select.style.display == 'block')
        select.style.display = 'none';

}


function show(an, width, height, borderStyle) {
    //selected();
    hide(href);
    //hide("select");

    //selected();
    var href = an.href;
    // alert(href);
    var boxdiv = document.getElementById(href);
    //var selected = document.getElementById("select");
    //	if (boxdiv != null) {
    //		if (boxdiv.style.display == 'none') {
    //           // move_box(an,boxdiv);
    //			fixed_box(boxdiv, 800, 50);
    //			boxdiv.style.display = 'block';
    //		} else
    //			boxdiv.style.display = 'none';
    //
    //		return false;
    //	}

    boxdiv = document.createElement('div');

    boxdiv.setAttribute('id', href);
    boxdiv.style.display = 'block';
    boxdiv.style.position = 'absolute';
    boxdiv.style.width = width + 'px';
    boxdiv.style.height = height + 'px';
    boxdiv.style.border = borderStyle;
    boxdiv.style.backgroundColor = '#fff';

    var contents = document.createElement('iframe');
    contents.scrolling = 'no';
    // contents.autoHideEnabled='yes'
    contents.frameBorder = '0';
    contents.style.width = width + 'px';
    contents.style.height = height + 'px';


    contents.src = href;

    boxdiv.appendChild(contents);
    document.body.appendChild(boxdiv);
    fixed_box(boxdiv, 800, 50);


    //selected = document.createElement('div');

    //selected.setAttribute('id', "select");

    //
    //	selected.style.display = 'block';
    //	selected.style.position = 'absolute';
    //
    //	selected.style.width = 20 + 'px';
    //	selected.style.height = 10 + 'px';
    //	selected.style.backgroundColor = '#fff';
    //	selected.innerHTML = "<html><table><tr><td  class=\"webSiteBodyGreenText\"><img src=\"/glycomics/common/images/red_arrow_up.gif\" width=\"12\" height=\"16\"/></td></td></table><html>";
    //	document.body.appendChild(selected);
    //move_box(an, selected, -10, 5);

    return false;
}

function SwitchMenu(obj) {
    if (document.getElementById) { //DynamicDrive.com change
        document.write('<style \n')
        document.write('.submenu{display: none;}\n')
        document.write('</style>\n')
    }

    if (document.getElementById) {
        var el = document.getElementById(obj);
        var ar = document.getElementById("masterdiv").getElementsByTagName("span");
        //DynamicDrive.com change
        if (el.style.display != "block") { //DynamicDrive.com change

            if (ar.className == "submenu") //DynamicDrive.com change
                ar.style.display = "none";

            el.style.display = "block";
        } else {
            el.style.display = "none";
        }
    }
}


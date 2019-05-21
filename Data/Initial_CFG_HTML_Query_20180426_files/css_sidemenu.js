
unhighlighted_bullet= new Image();
unhighlighted_bullet.src = "/static/consortium/images/samp2_bullet.gif";

highlighted_bullet = new Image();
highlighted_bullet.src = "/static/consortium/images/samp2_bullet_hl.gif";

plus = new Image();
plus.src = "/static/consortium/images/sub_plus.gif";

minus = new Image();
minus.src = "/static/consortium/images/sub_minus.gif";

function highlightBullet(bulletID) {
    bulletID.src=highlighted_bullet.src
}

function unhighlightBullet(bulletID) {
    bulletID.src=unhighlighted_bullet.src
}

function plusToMinus(imgID) {
    imgID.src=minus.src
}

function minusToPlus(imgID) {
    imgID.src=plus.src
}
/****************************************************
*	        DOM Image rollover:

*		by Chris Poole
*		http://chrispoole.com
*               Script featured on http://www.dynamicdrive.com
*		Keep this notice intact to use it :-)
****************************************************/
function init() {
  if (!document.getElementById) return
  var imgOriginSrc;
  var imgTemp = new Array();
  var imgarr = document.getElementsByTagName('img');
  for (var i = 0; i < imgarr.length; i++) {
    if (imgarr[i].getAttribute('hsrc')) {
        imgTemp[i] = new Image();
        imgTemp[i].src = imgarr[i].getAttribute('hsrc');
        imgarr[i].onmouseover = function() {
            imgOriginSrc = this.getAttribute('src');
            this.setAttribute('src',this.getAttribute('hsrc'))
        }
        imgarr[i].onmouseout = function() {
            this.setAttribute('src',imgOriginSrc)
        }
    }
  }
}
onload=init;
/***********************************************
* Switch Menu script- by Martial B of http://getElementById.com/
* Modified by Dynamic Drive for format & NS4/IE4 compatibility
* Visit http://www.dynamicdrive.com/ for full source code
***********************************************/
//if (document.getElementById){ //DynamicDrive.com change
  // document.write('<style type="text/css">')
   //document.write('.submenu{display: none;}')
   //document.write('</style>')
//}
function SwitchMenu(obj){
    if(document.getElementById){
        var el = document.getElementById(obj);
        var ar = document.getElementById("masterdiv").getElementsByTagName("span"); //DynamicDrive.com change
        if(el.style.display != "block"){ //DynamicDrive.com change
            for (var i=0; i<ar.length; i++){
                if (ar[i].className=="submenu") //DynamicDrive.com change
                    ar[i].style.display = "none";
            }
           el.style.display = "block";
       }else{
           el.style.display = "none";
       }
   }
}
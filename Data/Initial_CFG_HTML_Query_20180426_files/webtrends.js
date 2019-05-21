// WebTrends SmartSource Data Collector Tag
// Version: 8.6.0     
// Tag Builder Version: 2.1.0
// Created: 04/17/2009 4:49:40
// Updated: 04/06/2010

function WebTrends(){
	var that=this;
	// begin: user modifiable
	this.dcsid="dcsc2bh0q10000cprjxnf6sz5_6h4b";
	this.domain="statse.webtrendslive.com";
	this.timezone=-5;
	this.fpcdom="";
	this.onsitedoms="www.nature.com, blog-msb.embo.org, blogs.nature.com,cmckb.cellmigration.org, connotea.org, functionalglycomics.com, go.nature.com, guide.labanimal.com, kb.psi-structuralgenomics.org, lipidmaps.org, nature.com, natureasia.com, naturechina.com.cn, network.nature.com, neuroscience-gateway.org, pid.nci.nih.gov, precedings.nature.com, scintilla.nature.com, secure.nature.com, signaling-gateway.org, spotlight.nature.com, npg.nature.com, press.nature.com , palgrave-connect.com, www.cellmigration.org, www.connotea.org, www.econolog.net, www.functionalglycomics.org, www.labanimal.com, www.lipidmaps.org, www.natureasia.com, www.neuroscience-gateway.org, www.palgrave-journals.com, www.rikenresearch.riken.jp, www.scitable.com, www.signaling-gateway.org, cmkb.cellmigration.org, nrgwiki.nature.com, www.iaor-palgrave.com, www.natureprotocols.com, media.nature.com";
	this.downloadtypes="xls,doc,pdf,ris,txt,csv,zip,mp3,mp4,mov,AVI,bin,bmp,BMP,cdx,cml,CSS,css,dat,db,exe,Exe,GIF,gif,gz,HTM,jar,java,jpg,mol,pkg,png,ppt,psd,rss,sgml,tab,TAB,tgz,ttf,TXT,utf,xml,rdf";
	this.metanames="Access, citation_title, citation_journal_title, citation_date,citation_volume,citation_issue,citation_doi";
	this.trackevents=true;
	this.enabled=true;
	this.i18n=false;
	this.fpc="WT_FPC";
	// end: user modifiable
	this.DCS={};
	this.WT={};
	this.DCSext={};
	this.images=[];
	this.index=0;
	this.qp=[];
	this.exre=(function(){return(window.RegExp?new RegExp("dcs(uri)|(ref)|(aut)|(met)|(sta)|(sip)|(pro)|(byt)|(dat)|(p3p)|(cfg)|(redirect)|(cip)","i"):"");})();
	this.re=(function(){return(window.RegExp?(that.i18n?{"%25":/\%/g}:{"%09":/\t/g,"%20":/ /g,"%23":/\#/g,"%26":/\&/g,"%2B":/\+/g,"%3F":/\?/g,"%5C":/\\/g,"%22":/\"/g,"%7F":/\x7F/g,"%A0":/\xA0/g}):"");})();
}
WebTrends.prototype.dcsGetId=function(){
	if (this.enabled&&(document.cookie.indexOf(this.fpc+"=")==-1)&&(document.cookie.indexOf("WTLOPTOUT=")==-1)){
		document.write("<scr"+"ipt type='text/javascript' src='"+"http"+(window.location.protocol.indexOf('https:')==0?'s':'')+"://"+this.domain+"/"+this.dcsid+"/wtid.js"+"'><\/scr"+"ipt>");
	}
}
WebTrends.prototype.dcsGetCookie=function(name){
	var cookies=document.cookie.split("; ");
	var cmatch=[];
	var idx=0;
	var i=0;
	var namelen=name.length;
	var clen=cookies.length;
	for (i=0;i<clen;i++){
		var c=cookies[i];
		if ((c.substring(0,namelen+1))==(name+"=")){
			cmatch[idx++]=c;
		}
	}
	var cmatchCount=cmatch.length;
	if (cmatchCount>0){
		idx=0;
		if ((cmatchCount>1)&&(name==this.fpc)){
			var dLatest=new Date(0);
			for (i=0;i<cmatchCount;i++){
				var lv=parseInt(this.dcsGetCrumb(cmatch[i],"lv"));
				var dLst=new Date(lv);
				if (dLst>dLatest){
					dLatest.setTime(dLst.getTime());
					idx=i;
				}
			}
		}
		return unescape(cmatch[idx].substring(namelen+1));
	}
	else{
		return null;
	}
}
WebTrends.prototype.dcsGetCrumb=function(cval,crumb,sep){
	var aCookie=cval.split(sep||":");
	for (var i=0;i<aCookie.length;i++){
		var aCrumb=aCookie[i].split("=");
		if (crumb==aCrumb[0]){
			return aCrumb[1];
		}
	}
	return null;
}
WebTrends.prototype.dcsGetIdCrumb=function(cval,crumb){
	var id=cval.substring(0,cval.indexOf(":lv="));
	var aCrumb=id.split("=");
	for (var i=0;i<aCrumb.length;i++){
		if (crumb==aCrumb[0]){
			return aCrumb[1];
		}
	}
	return null;
}
WebTrends.prototype.dcsIsFpcSet=function(name,id,lv,ss){
	var c=this.dcsGetCookie(name);
	if (c){
		return ((id==this.dcsGetIdCrumb(c,"id"))&&(lv==this.dcsGetCrumb(c,"lv"))&&(ss==this.dcsGetCrumb(c,"ss")))?0:3;
	}
	return 2;
}
WebTrends.prototype.dcsFPC=function(){
	if (document.cookie.indexOf("WTLOPTOUT=")!=-1){
		return;
	}
	var WT=this.WT;
	var name=this.fpc;
	var dCur=new Date();
	var adj=(dCur.getTimezoneOffset()*60000)+(this.timezone*3600000);
	dCur.setTime(dCur.getTime()+adj);
	var dExp=new Date(dCur.getTime()+315360000000);
	var dSes=new Date(dCur.getTime());
	WT.co_f=WT.vtid=WT.vtvs=WT.vt_f=WT.vt_f_a=WT.vt_f_s=WT.vt_f_d=WT.vt_f_tlh=WT.vt_f_tlv="";
	if (document.cookie.indexOf(name+"=")==-1){
		if ((typeof(gWtId)!="undefined")&&(gWtId!="")){
			WT.co_f=gWtId;
		}
		else if ((typeof(gTempWtId)!="undefined")&&(gTempWtId!="")){
			WT.co_f=gTempWtId;
			WT.vt_f="1";
		}
		else{
			WT.co_f="2";
			var curt=dCur.getTime().toString();
			for (var i=2;i<=(32-curt.length);i++){
				WT.co_f+=Math.floor(Math.random()*16.0).toString(16);
			}
			WT.co_f+=curt;
			WT.vt_f="1";
		}
		if (typeof(gWtAccountRollup)=="undefined"){
			WT.vt_f_a="1";
		}
		WT.vt_f_s=WT.vt_f_d="1";
		WT.vt_f_tlh=WT.vt_f_tlv="0";
	}
	else{
		var c=this.dcsGetCookie(name);
		var id=this.dcsGetIdCrumb(c,"id");
		var lv=parseInt(this.dcsGetCrumb(c,"lv"));
		var ss=parseInt(this.dcsGetCrumb(c,"ss"));
		if ((id==null)||(id=="null")||isNaN(lv)||isNaN(ss)){
			return;
		}
		WT.co_f=id;
		var dLst=new Date(lv);
		WT.vt_f_tlh=Math.floor((dLst.getTime()-adj)/1000);
		dSes.setTime(ss);
		if ((dCur.getTime()>(dLst.getTime()+1800000))||(dCur.getTime()>(dSes.getTime()+28800000))){
			WT.vt_f_tlv=Math.floor((dSes.getTime()-adj)/1000);
			dSes.setTime(dCur.getTime());
			WT.vt_f_s="1";
		}
		if ((dCur.getDay()!=dLst.getDay())||(dCur.getMonth()!=dLst.getMonth())||(dCur.getYear()!=dLst.getYear())){
			WT.vt_f_d="1";
		}
	}
	WT.co_f=escape(WT.co_f);
	WT.vtid=(typeof(this.vtid)=="undefined")?WT.co_f:(this.vtid||"");
	WT.vtvs=(dSes.getTime()-adj).toString();
	var expiry="; expires="+dExp.toGMTString();
	var cur=dCur.getTime().toString();
	var ses=dSes.getTime().toString();
	document.cookie=name+"="+"id="+WT.co_f+":lv="+cur+":ss="+ses+expiry+"; path=/"+(((this.fpcdom!=""))?("; domain="+this.fpcdom):(""));
	var rc=this.dcsIsFpcSet(name,WT.co_f,cur,ses);
	if (rc!=0){
		WT.co_f=WT.vtvs=WT.vt_f_s=WT.vt_f_d=WT.vt_f_tlh=WT.vt_f_tlv="";
		if (typeof(this.vtid)=="undefined"){
			WT.vtid="";
		}
		WT.vt_f=WT.vt_f_a=rc;
    }
}
// Code section for Assign your query parameters to WebTrends query parameters. You will also need to configure the tag to recognize your query parameters.
WebTrends.prototype.dcsQP=function(N){
	if (typeof(N)=="undefined"){
		return "";
	}
	var qry=location.search.substring(1);
	if (qry!=""){
		var pairs=qry.split("&");
		for (var i=0;i<pairs.length;i++){
			var pos=pairs[i].indexOf("=");
			if (pos!=-1){
				if (pairs[i].substring(0,pos)==N){
					this.qp[this.qp.length]=(i==0?"":"&")+pairs[i];
					return pairs[i].substring(pos+1);
				}
			}
		}
	}
	return "";
}
WebTrends.prototype.dcsIsOnsite=function(host){
	if (host.length>0){
	    host=host.toLowerCase();
	    if (host==this.getHost()){
		    return true;
	    }
	    if (this.isonsiteHelper(host))
	    {
	    	return true;
	    }
	    if (typeof(this.onsitedoms.test)=="function"){
		    return this.onsitedoms.test(host);
	    }
	    else if (this.onsitedoms.length>0){
		    var doms=this.dcsSplit(this.onsitedoms);
		    var len=doms.length;
		    for (var i=0;i<len;i++){
			    if (host==doms[i]){
			        return true;
			    }
		    }
	    }
	}
	return false;
}
WebTrends.prototype.dcsTypeMatch=function(pth, typelist){
	var type=pth.substring(pth.lastIndexOf(".")+1,pth.length);
	var types=this.dcsSplit(typelist);
	var tlen=types.length;	
	for (var i=0;i<tlen;i++){
		if (type==types[i]){
			return true;
		}
	}
	return false;
}
WebTrends.prototype.dcsEvt=function(evt,tag){
	var e=evt.target||evt.srcElement;
	while (e.tagName&&(e.tagName!=tag)){
		e=e.parentElement||e.parentNode;
	}
	return e;
}
WebTrends.prototype.dcsNavigation=function(evt){
	return "";
}
WebTrends.prototype.dcsBind=function(event,func){
	if ((typeof(func)=="function")&&document.body){
		if (document.body.addEventListener){
			document.body.addEventListener(event, func.wtbind(this), true);
		}
		else if(document.body.attachEvent){
			document.body.attachEvent("on"+event, func.wtbind(this));
		}
	}
}
WebTrends.prototype.dcsET=function(){
	var e=(navigator.appVersion.indexOf("MSIE")!=-1)?"click":"click";
	this.dcsBind(e,this.dcsDownload);
	this.dcsBind(e,this.dcsOffsite);
}
WebTrends.prototype.dcsMultiTrack=function(){
	var args=dcsMultiTrack.arguments?dcsMultiTrack.arguments:arguments;
	if (args.length%2==0){
		this.dcsSetProps(args);
		var dCurrent=new Date();
		this.DCS.dcsdat=dCurrent.getTime();
		this.dcsFPC();
		this.dcsTag();
	}

}
WebTrends.prototype.dcsCleanUp=function(){
	this.DCS={};
	this.WT={};
	this.DCSext={};
	if (arguments.length%2==0){
		this.dcsSetProps(arguments);
	}
}
WebTrends.prototype.dcsSetProps=function(args){
	for (var i=0;i<args.length;i+=2){
		if (args[i].indexOf('WT.')==0){
			this.WT[args[i].substring(3)]=args[i+1];
		}
		else if (args[i].indexOf('DCS.')==0){
			this.DCS[args[i].substring(4)]=args[i+1];
		}
		else if (args[i].indexOf('DCSext.')==0){
			this.DCSext[args[i].substring(7)]=args[i+1];
		}
	}
}
WebTrends.prototype.dcsSplit=function(list){
	var items=list.toLowerCase().split(",");
	var len=items.length;
	for (var i=0;i<len;i++){
		items[i]=items[i].replace(/^\s*/,"").replace(/\s*$/,"");
	}
	return items;
}
// Code section for Track clicks to download links.
WebTrends.prototype.dcsDownload=function(evt){
	evt=evt||(window.event||"");
	if (evt&&((typeof(evt.which)!="number")||(evt.which==1))){
		var e=this.dcsEvt(evt,"A");
		if (e.href){
		    var hn=e.hostname?(e.hostname.split(":")[0]):"";
		    if (this.dcsIsOnsite(hn)&&this.dcsTypeMatch(e.pathname,this.downloadtypes)){
		        var qry=e.search?e.search.substring(e.search.indexOf("?")+1,e.search.length):"";
		        var pth=e.pathname?((e.pathname.indexOf("/")!=0)?"/"+e.pathname:e.pathname):"/";
		        var ttl="";
		        var text=document.all?e.innerText:e.text;
		        var img=this.dcsEvt(evt,"IMG");
		        if (img.alt){
			        ttl=img.alt;
		        }
		        else if (text){
			        ttl=text;
		        }
		        else if (e.innerHTML){
			        ttl=e.innerHTML;
		        }
		        this.dcsMultiTrack("DCS.dcssip",hn,"DCS.dcsuri",pth,"DCS.dcsqry",e.search||"","WT.ti","Download:"+ttl,"WT.dl","20","WT.ndl","0","WT.nv",this.dcsNavigation(evt));
		        this.DCS.dcssip=this.DCS.dcsuri=this.DCS.dcsqry=this.WT.ti=this.WT.dl=this.WT.nv=this.WT.ndl="";
		    }
		}
	}
}
// Code section for Track clicks to links leading offsite.
WebTrends.prototype.dcsOffsite=function(evt){
	evt=evt||(window.event||"");
	if (evt&&((typeof(evt.which)!="number")||(evt.which==1))){
		var e=this.dcsEvt(evt,"A");
		if (e.href){
		    var hn=e.hostname?(e.hostname.split(":")[0]):"";
		    var pr=e.protocol||"";
		    if ((hn.length>0)&&(pr.indexOf("http")==0)&&!this.dcsIsOnsite(hn)){
			    var qry=e.search?e.search.substring(e.search.indexOf("?")+1,e.search.length):"";
			    var pth=e.pathname?((e.pathname.indexOf("/")!=0)?"/"+e.pathname:e.pathname):"/";
			    var trim=true;
			    this.dcsMultiTrack("DCS.dcssip",hn,"DCS.dcsuri",pth,"DCS.dcsqry",trim?"":qry,"WT.ti","Offsite:"+hn+pth+qry,"WT.dl","24","WT.ndl","24","WT.nv",this.dcsNavigation(evt));
			    this.DCS.dcssip=this.DCS.dcsuri=this.DCS.dcsqry=this.WT.ti=this.WT.dl=this.WT.nv=this.WT.ndl="";
		    }
		}
	}
}

// Code section for Enable custom meta tag capture.
WebTrends.prototype.dcsMetaCap=function(){
	var elems;
	if (document.all){
		elems=document.all.tags("meta");
	}
	else if (document.documentElement){
		elems=document.getElementsByTagName("meta");
	}
	if (typeof(elems)!="undefined"){
		var names=this.dcsSplit(this.metanames);
		var nlen=names.length;
		var elen=elems.length;
		for (var i=0;i<elen;i++){
			var name=elems.item(i).name.toLowerCase();
			var content=elems.item(i).content;
			if ((name.length>0)&&(content.length>0)){
				for(var j=0;j<nlen;j++){
					if (name==names[j]){
						this.DCSext["meta_"+names[j]]=content;
						break;
					}
				}
			}
		}
	}
}

WebTrends.prototype.dcsAdv=function(){
	if (this.trackevents&&(typeof(this.dcsET)=="function")){
		if (window.addEventListener){
			window.addEventListener("load",this.dcsET.wtbind(this),false);
		}
		else if (window.attachEvent){
			window.attachEvent("onload",this.dcsET.wtbind(this));
		}
	}
	this.dcsFPC();
	this.dcsMetaCap();
}
WebTrends.prototype.dcsVar=function(){
	var dCurrent=new Date();
	var WT=this.WT;
	var DCS=this.DCS;
	WT.tz=parseInt(dCurrent.getTimezoneOffset()/60*-1)||"0";
	WT.bh=dCurrent.getHours()||"0";
	WT.ul=navigator.appName=="Netscape"?navigator.language:navigator.userLanguage;
	if (typeof(screen)=="object"){
		WT.cd=navigator.appName=="Netscape"?screen.pixelDepth:screen.colorDepth;
		WT.sr=screen.width+"x"+screen.height;
	}
	if (typeof(navigator.javaEnabled())=="boolean"){
		WT.jo=navigator.javaEnabled()?"Yes":"No";
	}
	if (document.title){
		if (window.RegExp){
			var tire=new RegExp("^"+window.location.protocol+"//"+this.getHost()+"\\s-\\s");
			WT.ti=document.title.replace(tire,"");
		}
		else{
			WT.ti=document.title;
		}
	}
	WT.js="Yes";
	WT.jv=(function(){
		var agt=navigator.userAgent.toLowerCase();
		var major=parseInt(navigator.appVersion);
		var mac=(agt.indexOf("mac")!=-1);
		var ff=(agt.indexOf("firefox")!=-1);
		var ff0=(agt.indexOf("firefox/0.")!=-1);
		var ff10=(agt.indexOf("firefox/1.0")!=-1);
		var ff15=(agt.indexOf("firefox/1.5")!=-1);
		var ff20=(agt.indexOf("firefox/2.0")!=-1);
		var ff3up=(ff&&!ff0&&!ff10&!ff15&!ff20);
		var nn=(!ff&&(agt.indexOf("mozilla")!=-1)&&(agt.indexOf("compatible")==-1));
		var nn4=(nn&&(major==4));
		var nn6up=(nn&&(major>=5));
		var ie=((agt.indexOf("msie")!=-1)&&(agt.indexOf("opera")==-1));
		var ie4=(ie&&(major==4)&&(agt.indexOf("msie 4")!=-1));
		var ie5up=(ie&&!ie4);
		var op=(agt.indexOf("opera")!=-1);
		var op5=(agt.indexOf("opera 5")!=-1||agt.indexOf("opera/5")!=-1);
		var op6=(agt.indexOf("opera 6")!=-1||agt.indexOf("opera/6")!=-1);
		var op7up=(op&&!op5&&!op6);
		var jv="1.1";
		if (ff3up){
			jv="1.8";
		}
		else if (ff20){
			jv="1.7";
		}
		else if (ff15){
			jv="1.6";
		}
		else if (ff0||ff10||nn6up||op7up){
			jv="1.5";
		}
		else if ((mac&&ie5up)||op6){
			jv="1.4";
		}
		else if (ie5up||nn4||op5){
			jv="1.3";
		}
		else if (ie4){
			jv="1.2";
		}
		return jv;
	})();
	WT.ct="unknown";
	if (document.body&&document.body.addBehavior){
		try{
			document.body.addBehavior("#default#clientCaps");
			WT.ct=document.body.connectionType||"unknown";
			document.body.addBehavior("#default#homePage");
			WT.hp=document.body.isHomePage(location.href)?"1":"0";
		}
		catch(e){
		}
	}
	if (document.all){
		WT.bs=document.body?document.body.offsetWidth+"x"+document.body.offsetHeight:"unknown";
	}
	else{
		WT.bs=window.innerWidth+"x"+window.innerHeight;
	}
	WT.fv=(function(){
		var i,flash;
		if (window.ActiveXObject){
			for(i=10;i>0;i--){
				try{
					flash=new ActiveXObject("ShockwaveFlash.ShockwaveFlash."+i);
					return i+".0";
				}
				catch(e){
				}
			}
		}
		else if (navigator.plugins&&navigator.plugins.length){
			for (i=0;i<navigator.plugins.length;i++){
				if (navigator.plugins[i].name.indexOf('Shockwave Flash')!=-1){
					return navigator.plugins[i].description.split(" ")[2];
				}
			}
		}
		return "Not enabled";
	})();
	WT.slv=(function(){
		var slv="Not enabled";
		try{     
			if (navigator.userAgent.indexOf('MSIE')!=-1){
				var sli = new ActiveXObject('AgControl.AgControl');
				if (sli){
					slv="Unknown";
				}
			}
			else if (navigator.plugins["Silverlight Plug-In"]){
				slv="Unknown";
			}
		}
		catch(e){
		}
		if (slv!="Not enabled"){
			var i,j,v;
			if ((typeof(Silverlight)=="object")&&(typeof(Silverlight.isInstalled)=="function")){
				for (j=9;j>=0;j--){
					for(i=3;i>0;i--){
						v=i+"."+j;
						if (Silverlight.isInstalled(v)){
							slv=v;
							break;
						}
					}
					if (slv==v){
						break;
					}
				}
			}
		}
		return slv;
	})();
	if (this.i18n){
		if (typeof(document.defaultCharset)=="string"){
			WT.le=document.defaultCharset;
		} 
		else if (typeof(document.characterSet)=="string"){
			WT.le=document.characterSet;
		}
		else{
			WT.le="unknown";
		}
	}
	WT.tv="8.6.0";
//	WT.sp="@@SPLITVALUE@@";
	WT.dl="0";
	WT.ndl="0";
	WT.ssl=(window.location.protocol.indexOf('https:')==0)?"1":"0";
	DCS.dcsdat=dCurrent.getTime();
	DCS.dcssip=this.getHost();
	DCS.dcsuri=window.location.pathname;
	WT.es=DCS.dcssip+DCS.dcsuri;
	if (window.location.search){
		DCS.dcsqry=window.location.search;
		if (this.qp.length>0){
			for (var i=0;i<this.qp.length;i++){
				var pos=DCS.dcsqry.indexOf(this.qp[i]);
				if (pos!=-1){
					var front=DCS.dcsqry.substring(0,pos);
					var end=DCS.dcsqry.substring(pos+this.qp[i].length,DCS.dcsqry.length);
					DCS.dcsqry=front+end;
				}
			}
		}
	}
	if ((window.document.referrer!="")&&(window.document.referrer!="-")){
		if (!(navigator.appName=="Microsoft Internet Explorer"&&parseInt(navigator.appVersion)<4)){
			DCS.dcsref=window.document.referrer;
		}
	}
}
WebTrends.prototype.dcsEscape=function(S, REL){
	if (REL!=""){
		S=S.toString();
		for (var R in REL){
 			if (REL[R] instanceof RegExp){
				S=S.replace(REL[R],R);
 			}
		}
		return S;
	}
	else{
		return escape(S);
	}
}
WebTrends.prototype.dcsA=function(N,V){
	if (this.i18n&&(this.exre!="")&&!this.exre.test(N)){
		if (N=="dcsqry"){
			var newV="";
			var params=V.substring(1).split("&");
			for (var i=0;i<params.length;i++){
				var pair=params[i];
				var pos=pair.indexOf("=");
				if (pos!=-1){
					var key=pair.substring(0,pos);
					var val=pair.substring(pos+1);
					if (i!=0){
						newV+="&";
					}
					newV+=key+"="+this.dcsEncode(val);
				}
			}
			V=V.substring(0,1)+newV;
		}
		else{
			V=this.dcsEncode(V);
		}
	}
	return "&"+N+"="+this.dcsEscape(V, this.re);
}
WebTrends.prototype.dcsEncode=function(S){
	return (typeof(encodeURIComponent)=="function")?encodeURIComponent(S):escape(S);
}
WebTrends.prototype.dcsCreateImage=function(dcsSrc){
	if (document.images){
		this.images[this.index]=new Image();
		this.images[this.index].src=dcsSrc;
		this.index++;
	}
	else{
		document.write('<IMG ALT="" BORDER="0" NAME="DCSIMG" WIDTH="1" HEIGHT="1" SRC="'+dcsSrc+'">');
	}
}
WebTrends.prototype.dcsMeta=function(){
	var elems;
	if (document.all){
		elems=document.all.tags("meta");
	}
	else if (document.documentElement){
		elems=document.getElementsByTagName("meta");
	}
	if (typeof(elems)!="undefined"){
		var length=elems.length;
		for (var i=0;i<length;i++){
			var name=elems.item(i).name;
			var content=elems.item(i).content;
			var equiv=elems.item(i).httpEquiv;
			if (name.length>0){
				if (name.toUpperCase().indexOf("WT.")==0){
					this.WT[name.substring(3)]=content;
				}
				else if (name.toUpperCase().indexOf("DCSEXT.")==0){
					this.DCSext[name.substring(7)]=content;
				}
				else if (name.toUpperCase().indexOf("DCS.")==0){
					this.DCS[name.substring(4)]=content;
				}
			}
		}
	}
}
WebTrends.prototype.dcsTag=function(){
	if (document.cookie.indexOf("WTLOPTOUT=")!=-1){
		return;
	}
	var WT=this.WT;
	var DCS=this.DCS;
	var DCSext=this.DCSext;
	var i18n=this.i18n;
	var P="http"+(window.location.protocol.indexOf('https:')==0?'s':'')+"://"+this.domain+(this.dcsid==""?'':'/'+this.dcsid)+"/dcs.gif?";
	if (i18n){
		WT.dep="";
	}
	for (var N in DCS){
 		if (DCS[N]&&(typeof DCS[N]!="function")){
			P+=this.dcsA(N,DCS[N]);
		}
	}
	var keys=["co_f","vtid","vtvs","vt_f_tlv"];
	for (var i=0;i<keys.length;i++){
		var key=keys[i];
		if (WT[key]){
			P+=this.dcsA("WT."+key,WT[key]);
			delete WT[key];
		}
	}
	for (N in WT){
		if (WT[N]&&(typeof WT[N]!="function")){
			P+=this.dcsA("WT."+N,WT[N]);
		}
	}
	for (N in DCSext){
		if (DCSext[N]&&(typeof DCSext[N]!="function")){
			if (i18n){
				WT.dep=(WT.dep.length==0)?N:(WT.dep+";"+N);
			}
			P+=this.dcsA(N,DCSext[N]);
		}
	}
	if (i18n&&(WT.dep.length>0)){
		P+=this.dcsA("WT.dep",WT.dep);
	}
	if (P.length>2048&&navigator.userAgent.indexOf('MSIE')>=0){
		P=P.substring(0,2040)+"&WT.tu=1";
	}
	this.dcsCreateImage(P);
	this.WT.ad="";
}
WebTrends.prototype.dcsDebug=function(){
	var t=this;
	var i=t.images[0].src;
	var q=i.indexOf("?");
	var r=i.substring(0,q).split("/");
	var m="<b>Protocol</b><br><code>"+r[0]+"<br></code>";
	m+="<b>Domain</b><br><code>"+r[2]+"<br></code>";
	m+="<b>Path</b><br><code>/"+r[3]+"/"+r[4]+"<br></code>";
	m+="<b>Query Params</b><code>"+i.substring(q+1).replace(/\&/g,"<br>")+"</code>";
	m+="<br><b>Cookies</b><br><code>"+document.cookie.replace(/\;/g,"<br>")+"</code>";
	if (t.w&&!t.w.closed){
		t.w.close();

	}
	t.w=window.open("","dcsDebug","width=500,height=650,scrollbars=yes,resizable=yes");
	t.w.document.write(m);
	t.w.focus();
}

WebTrends.prototype.isonsiteHelper=function(host)
{
    var osda=this.onsitedoms.split(",");
    var len=osda.length;
    for (x=0; x<len; x++)
    {
    	var o=osda[x].trim();
    	if (host.indexOf(o.toLowerCase())==0)
    	{
    		return true;
    	}
    }
    return false;	
}
WebTrends.prototype.dcsCollect=function(){
    if (this.enabled){
        this.dcsVar();
        this.dcsMeta();
        this.dcsAdv();
        this.dcsTag();
    }
}

//custom function to get window host
WebTrends.prototype.getHost=function(){
    host = window.location.hostname.toLowerCase();
    //if host starts with www.nature.com but also includes gatway domain (is longer than 14 char)
    if ((host.substr(0,14)=="www.nature.com")&&(host.length>14)) {
        _tag.DCSext.gateway=host.substr(14,host.length);
        //return trimed version of host name
        return host.substr(0,14);
    } else {
        //else return full host
        return host;
    }
}

function dcsMultiTrack(){
	if (typeof(_tag)!="undefined"){
		return(_tag.dcsMultiTrack());
	}
}

function dcsDebug(){
	if (typeof(_tag)!="undefined"){
		return(_tag.dcsDebug());
	}
}

Function.prototype.wtbind = function(obj){
	var method=this;
	var temp=function(){
		return method.apply(obj,arguments);
	};
	return temp;
}

String.prototype.trim = function () {
    return this.replace(/^\s*/, "").replace(/\s*$/, "");
}

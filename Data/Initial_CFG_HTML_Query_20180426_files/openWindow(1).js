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
function openWindow(page) {
    if (!popupWin || popupWin.closed) {
        popupWin = window.open(page, uniqueWindow, "width=1200,height=600,resizable,scrollbars");

    }
    else {
        popupWin = window.open(page, uniqueWindow, "width=1200,height=600,resizable,scrollbars");

        popupWin.focus();
    }
}

var Pic = new Array

Pic[0] = '/glycomics/common/images/presentation/expt_active.jpg'
Pic[1] = '/glycomics/common/images/presentation/expt_inactive.jpg'
Pic[2] = '/glycomics/common/images/presentation/summary_active.jpg'
Pic[3] = '/glycomics/common/images/presentation/summary_inactive.jpg'
Pic[4] = '/glycomics/common/images/presentation/data_active.jpg'
Pic[5] = '/glycomics/common/images/presentation/data_inactive.jpg'
Pic[6] = '/glycomics/common/images/presentation/summary_active.jpg'
Pic[7] = '/glycomics/common/images/presentation/summary_inactive.jpg'
Pic[8] = '/glycomics/common/images/presentation/array_active.jpg'
Pic[9] = '/glycomics/common/images/presentation/array_inactive.jpg'
Pic[10] = '/glycomics/common/images/presentation/excel_active.jpg'
Pic[11] = '/glycomics/common/images/presentation/excel_inactive.jpg'
Pic[12] = '/glycomics/common/images/presentation/data_active.jpg'
Pic[13] = '/glycomics/common/images/presentation/data_inactive.jpg'
Pic[14] = '/glycomics/common/images/presentation/search_active.gif'
Pic[15] = '/glycomics/common/images/presentation/search_inactive.gif'
Pic[16] = '/glycomics/common/images/presentation/data_active.jpg'
Pic[17] = '/glycomics/common/images/presentation/data_inactive.jpg'
Pic[18] = '/glycomics/common/images/presentation/rma_active.gif'
Pic[19] = '/glycomics/common/images/presentation/rma_inactive.gif'
Pic[20] = '/glycomics/common/images/presentation/download_active.gif'
Pic[21] = '/glycomics/common/images/presentation/download_inactive.gif'
Pic[22] = '/glycomics/common/images/presentation/link_active.gif'
Pic[23] = '/glycomics/common/images/presentation/link_inactive.gif'
Pic[24] = '/glycomics/common/images/presentation/rma_active.gif'
Pic[25] = '/glycomics/common/images/presentation/rma_inactive.gif'
Pic[26] = '/glycomics/common/images/presentation/link_active.gif'
Pic[27] = '/glycomics/common/images/presentation/link_inactive.gif'
Pic[28] = '/glycomics/common/images/presentation/cdf_active.gif'
Pic[29] = '/glycomics/common/images/presentation/cdf_inactive.gif'
Pic[30] = '/glycomics/common/images/presentation/desc_active.gif'
Pic[31] = '/glycomics/common/images/presentation/desc_inactive.gif'


var p = Pic.length
var preLoadImg = new Array()
var i = 0;
for (i = 0; i < p; i++) {
    preLoadImg[i] = new Image()
    preLoadImg[i].src = Pic[i]
}


function switchPresentationImage(whichImage, imageNumber) {
    document.images[whichImage].src = preLoadImg[imageNumber].src
}


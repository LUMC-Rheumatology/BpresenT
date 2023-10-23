function clearInput(form){
    var inputelements = form.getElementsByTagName("input");
    for (var ii=0; ii < inputelements.length; ii++) {
      if (inputelements[ii].type == "text") {
        inputelements[ii].value = "";
      }
    }
    var selectelements = form.getElementsByTagName("select");
    for (var jj=0; jj < selectelements.length; jj++) {
        selectelements[jj].value = null;
    }
}
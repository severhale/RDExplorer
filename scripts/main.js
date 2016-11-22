$("nav li").click(function (event) {
    event.preventDefault();
    var sectionId = $(this).attr('data-section');
    if (sectionId != 'simulations') {
        Processing.getInstanceById('sketch').noLoop();
    } else {
        Processing.getInstanceById('sketch').loop();
    }
    $(".active").removeClass('active');
    $("#" + sectionId).addClass('active');

    $(".selected").removeClass("selected");
    $(this).addClass("selected");
});


// why won't this bs work
// it goddamn works in regular processing
//function updateActivatorRadius(r) {
//    console.log("Attempting to set activator radius");
//    console.log(r);
//    var pInst = Processing.getInstanceById('sketch');
//    pInst.setActivatorRadius(r);
//}
//
//function updateInhibitorRadius(r) {
//    console.log("Attempting to set inhibitor radius");
//    console.log(r);
//    var pInst = Processing.getInstanceById('sketch');
//    pInst.setInhibitorRadius(r);
//}
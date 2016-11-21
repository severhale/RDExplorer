$("nav li").click(function (event) {
    event.preventDefault();
    var sectionId = $(this).attr('data-section');
    $(".active").removeClass('active');
    $("#" + sectionId).addClass('active');

    $(".selected").removeClass("selected");
    $(this).addClass("selected");
});
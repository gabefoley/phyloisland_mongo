// jQuery(document).ready(function () {
//
//     $('#submit_hit').on('click', function (e) {
//         e.preventDefault();
//         console.log('clicked submit hit')
//
//     });
//
//     $('#delete_hit').on('click', function (e) {
//         e.preventDefault();
//         console.log('clicked delete hit')
//
//         var items = "{{ items |safe }}"
//         hits = "gabe"
//
//         $.ajax({
//             url: Flask.url_for("delete_tag"),
//             data: JSON.stringify({
//                 'items' : items, 'hits' : hits
//             }),
//             contentType: 'application/json;charset=UTF-8',
//             type: 'POST',
//             error: function(error) {
//                 console.log(error);
//             }
//         });
//
//     });
//
// });
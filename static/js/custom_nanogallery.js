var baseURL = 'https://nanogallery2.nanostudio.org/samples/';
console.log('got here')
// add one image

                 jQuery(document).ready(function () {

                     jQuery('#btn_add').on('click', function () {
                         var ngy2data = $("#my_nanogallery2").nanogallery2('data');
                         var instance = $("#my_nanogallery2").nanogallery2('instance');

                         // create the new item
                         var ID = ngy2data.items.length + 1;
                         var albumID = '0';
                         var newItem = NGY2Item.New(instance, 'new berlin ' + ID, '', ID, albumID, 'image', '');

                         // define thumbnail -> image displayed in the gallery
                         newItem.thumbSet(baseURL + 'berlin2t.jpg', 135, 204); // w,h

                         // define URL of the media to display in lightbox
                         newItem.setMediaURL(baseURL + 'berlin2.jpg', 'img');

                         // currently displayed album
                         if (ngy2data.items[ngy2data.gallery.albumIdx].GetID() == 0) {

                             // add new item to current Gallery Object Model (GOM)
                             newItem.addToGOM();

                             // refresh the display (-> only once if you add multiple items)
                             $("#my_nanogallery2").nanogallery2('resize');
                         }
                     });


// add 2 images to existing album
                     jQuery('#btn_add_album').on('click', function () {
                         var ngy2data = $("#my_nanogallery2").nanogallery2('data');
                         var instance = $("#my_nanogallery2").nanogallery2('instance');

                         var ID = ngy2data.items.length + 1001;
                         var albumID = '1000';   // ID of the album which will contain the 2 new images

                         // ADD FIRST IMAGE
                         var newItem1 = NGY2Item.New(instance, 'new image ' + ID, '', ID, albumID, 'image', '');
                         newItem1.thumbSet(baseURL + 'berlin1t.jpg', 320, 212);
                         newItem1.setMediaURL(baseURL + 'berlin1.jpg', 'img');

                         // ADD SECOND IMAGE
                         ID = ngy2data.items.length + 1001;
                         var newItem2 = NGY2Item.New(instance, 'new image ' + ID, '', ID, albumID, 'image', '');
                         newItem2.thumbSet(baseURL + 'berlin1t.jpg', 320, 212);
                         newItem2.setMediaURL(baseURL + 'berlin1.jpg', 'img');

                         // refresh the display if the album is displayed
                         if (ngy2data.items[ngy2data.gallery.albumIdx].GetID() == albumID) {
                             newItem1.addToGOM();
                             newItem2.addToGOM();
                             $("#my_nanogallery2").nanogallery2('resize');
                         }
                         else {
                             alert('images added');
                         }
                     });


// delete selection
                     jQuery('#btn_del').on('click', function () {
                         console.log('was clicked')
                         var ngy2data = $("#my_nanogallery2").nanogallery2('data');
                         ngy2data.items.forEach(function (item) {
                             if (item.selected) {
                                 item.delete();
                             }
                         });
                         $("#my_nanogallery2").nanogallery2('resize');
                     });

// switch selection mode on/off
                     jQuery('#btn_select_mode').on('click', function () {
                         var b = !$('#my_nanogallery2').nanogallery2('option', 'thumbnailSelectable');
                         $('#my_nanogallery2').nanogallery2('option', 'thumbnailSelectable', b);
                     });

                 });
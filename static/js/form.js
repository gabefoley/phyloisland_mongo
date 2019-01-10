
    $(document).ready(function () {
        //attach the click event to all <img> tag having class 'popos_img'
        //here no need of any for loop
        $('.upload_submit').click(function () {
            $("#loading").show();
        });

        // Show / hide different parts of the form based on what type of file we're working with
        $('#type').change(function () {
            if ($(this).val() === 'species') {
                $("input[type='checkbox']").prop('disabled', false).prop('checked', true);
                $('.genome-search').css('color', 'black');
                $('#add_sequence').prop('disabled', true);
                $('#add_sequence').prop('checked', false);
                $("label[for='add_sequence']").css('color', '#ccc');
                $('#genome_type')[0].selectize.enable();

            }
            else if ($(this).val() === 'profile') {
                $("input[type='checkbox']").prop('disabled', true).prop('checked', false);
                $('.genome-search').css('color', '#ccc');
                $("label[for='add_sequence']").css('color', '#ccc');
                $('#genome_type')[0].selectize.disable();


            }

            else {
                $("input[type='checkbox']").prop('disabled', false).prop('checked', true);
                $('.genome-search').css('color', 'black');
                $("label[for='add_sequence']").css('color', 'black');
                $('#genome_type')[0].selectize.enable();



            }
        });

        // Drag and drop multi-selectable

        $('#input-draggable').selectize({
            plugins: ['drag_drop'],
            delimiter: ',',
            persist: false,
            create: function (input) {
                return {
                    value: input,
                    text: input
                };
            }
        });

        $('.input-sortable').selectize({
            plugins: ['drag_drop'],
            persist: false,
            create: true
        });

        $(function() {
          var $wrapper = $('#wrapper');



          // display scripts on the page
          $('script', $wrapper).each(function() {
            var code = this.text;
            if (code && code.length) {
              var lines = code.split('\n');
              var indent = null;

              for (var i = 0; i < lines.length; i++) {
                if (/^[  ]*$/.test(lines[i])) continue;
                if (!indent) {
                  var lineindent = lines[i].match(/^([  ]+)/);
                  if (!lineindent) break;
                  indent = lineindent[1];
                }
                lines[i] = lines[i].replace(new RegExp('^' + indent), '');
              }

              code = $.trim(lines.join('\n')).replace(/ /g, '    ');
              var $pre = $('<pre>').addClass('js').text(code);
              $pre.insertAfter(this);
            }
          });

          // show current input values
          $('select.selectized,input.selectized', $wrapper).each(function() {
            var $container = $('<div>').addClass('value').html('Current Value: ');
            var $value = $('<span>').appendTo($container);
            var $input = $(this);
            var update = function(e) { $value.text(JSON.stringify($input.val())); };

            $(this).on('change', update);
            update();

            $container.insertAfter($input);
          });
        });



        $('#genome_type').selectize({
            maxItems: null,
            valueField: 'id',
            labelField: 'title',
            searchField: 'title',
            plugins: ['drag_drop', 'remove_button'],
            create: false,
            highlight: true,


        });



        // Allow user to add another field
        $('.multi-field-wrapper').each(function() {
    var $wrapper = $('.multi-fields', this);
    $(".add-field", $(this)).click(function(e) {
        $('.multi-field:first-child', $wrapper).clone(true).appendTo($wrapper).find('input').val('').focus();
    });
    $('.multi-field .remove-field', $wrapper).click(function() {
        if ($('.multi-field', $wrapper).length > 1)
            $(this).parent('.multi-field').remove();
    });
});


    });







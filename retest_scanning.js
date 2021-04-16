var retest_data;
var bc_string = '';
var bcs_scanned = new Set();

function make_retest_scanning_table() {
  /**
   * Constructs the table from the main data
   */
  d3.select('#whoo_table_body').selectAll('.tr')
    .data(this.retest_data)
    .enter()
    .append('tr')
      .attr('class', 'retest_data_row');

  d3.selectAll('.retest_data_row').append('td').attr('class', 'Well_row_entry').html(d => d['Well']);
  d3.selectAll('.retest_data_row').append('td').attr('class', 'SampleName_row_entry').html(d => d['Sample_ID']);
  d3.selectAll('.retest_data_row').append('td').attr('class', 'Retest_row_entry').html(d => d['Result']);
  
  let pulled_entry = d3.selectAll('.retest_data_row').append('td').attr('class', 'Pulled_row_entry');
  d3.selectAll('.retest_data_row').append('td').attr('class', 'CT_row_entry').html(d => d['N1']);
  d3.selectAll('.retest_data_row').append('td').attr('class', 'CT_row_entry').html(d => d['RDRP']);
  d3.selectAll('.retest_data_row').append('td').attr('class', 'CT_row_entry').html(d => d['RNASEP']);
  pulled_entry.append('input')
    .attr('type', 'checkbox')
    .attr('class', 'pulled_checkbox')
    .on('change', function(e, d) {
      if (this.checked) {
        bcs_scanned.add(d['Sample_ID']);
      } else {
        bcs_scanned.delete(d['Sample_ID']);
      }
      d3.selectAll('.retest_data_row').classed('scanned_row', td => bcs_scanned.has(td['Sample_ID']));
    }); 

  d3.select('body').on('keydown', function(e) {
    if (e.key=='Clear') {
      bc_string = '';
    } else if (e.key=='Enter') {
      console.log('bc scanned:', bc_string);
      bcs_scanned.add(bc_string);
      d3.selectAll('.retest_data_row').classed('scanned_row', d => bcs_scanned.has(d['Sample_ID']));
      d3.selectAll('.pulled_checkbox').property('checked', function(d) {
        if (d['Sample_ID']==bc_string) {
          return true;
        } else {
          return this.checked;
        }
      });
    } else if (e.key != 'Shift') {
      bc_string += e.key;
    }
  });
}

function retest_scanning_setup() {
  let fr = new FileReader(); 
  fr.onload=function() { // this function runs once the file is loaded - all the action is in here
    d3.select('#inputfile').style('display', 'none'); // hide file input once a file is loaded
    d3.select('#plate_title').html('Retests file: ' + input_file.name); // display file name
    retest_data = d3.csvParse(fr.result);
    make_retest_scanning_table();
  }
  let input_file = document.getElementById('inputfile').files[0];
  fr.readAsText(document.getElementById('inputfile').files[0]);
}
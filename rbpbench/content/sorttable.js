

function isNumeric(value) {
    return /^-?\d+(\.\d+)?([eE][-+]?\d+)?$/.test(value);
}


function isNumericColumn(rows, columnIndex) {
    for (let i = 1; i < Math.min(3, rows.length); i++) {
        if (!isNumeric(rows[i].cells[columnIndex].textContent.trim())) {
            return false;
        }
    }
    return true;
}


function sortRows(rows, columnIndex, ascending, isNumericSort) {
    return rows.slice(1).sort(function(rowA, rowB) {
        var cellA = rowA.cells[columnIndex].textContent.trim();
        var cellB = rowB.cells[columnIndex].textContent.trim();

        if (isNumericSort) {
            cellA = parseFloat(cellA);
            cellB = parseFloat(cellB);
        }

        if (cellA < cellB) {
            return ascending ? -1 : 1;
        }
        if (cellA > cellB) {
            return ascending ? 1 : -1;
        }
        return 0;
    });
}


function rebuildTable(table, rows) {
    var newTbody = document.createElement("tbody");
    rows.forEach(function(row) {
        newTbody.appendChild(row);
    });

    var oldTbody = table.getElementsByTagName("tbody")[0];
    table.replaceChild(newTbody, oldTbody);
}


function makeSortable(table) {
    var headers = table.getElementsByTagName("th");
    var tableData = Array.from(table.getElementsByTagName("tr"));
    var sortDirections = new Array(headers.length).fill(true); // true for ascending

    for (var i = 0; i < headers.length; i++) {
        (function(index){
            headers[i].addEventListener('click', function() {
                var isNumericSort = isNumericColumn(tableData, index);
                var sortedRows = sortRows(tableData, index, sortDirections[index], isNumericSort);
                rebuildTable(table, sortedRows);
                sortDirections[index] = !sortDirections[index]; // toggle the direction
            });
        })(i);
    }
}


// Apply sorting to all tables in HTML document.
document.addEventListener('DOMContentLoaded', function() {
    var tables = document.getElementsByTagName('table');
    for (var i = 0; i < tables.length; i++) {
        makeSortable(tables[i]);
    }
});


















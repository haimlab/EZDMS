document.getElementById('fastaFile').addEventListener('change', function(e) {
    var file = this.files[0];
    var errorSpan = this.nextElementSibling;
    if (file) {
        var ext = file.name.split('.').pop().toLowerCase();
        if (ext !== "fasta" && ext !== "fa") {
            errorSpan.style.display = 'block';
            this.value = ''; // Clear the input
        } else {
            errorSpan.style.display = 'none';
        }
    }
});

document.getElementById('fastqFile').addEventListener('change', function(e) {
    var file = this.files[0];
    var errorSpan = this.nextElementSibling;
    if (file) {
        var ext = file.name.split('.').pop().toLowerCase();
        if (ext !== "fastq") {
            errorSpan.style.display = 'block';
            this.value = ''; // Clear the input
        } else {
            errorSpan.style.display = 'none';
        }
    }
});


document.addEventListener("DOMContentLoaded", function() {
    var helpSections = document.querySelectorAll('.container-help .help-box .help-section h2');
    helpSections.forEach(function(header) {
        header.addEventListener('click', function() {
            // Toggle visibility of nextElementSibling
            this.parentElement.classList.toggle('active');
            let content = this.nextElementSibling;
            if (content.style.display === "block") {
                content.style.display = "none";
            } else {
                content.style.display = "block";
            }
        });
    });
});



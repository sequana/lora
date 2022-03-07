/** add event listener to table row */
[].slice.call(document.querySelectorAll("#lora_summary tr"), 1).forEach(row => {
    const selectedClass = "is-selected";
    row.addEventListener("click", () => {
        /* get row informations */
        const ths = document.querySelectorAll("#lora_summary th", 1);
        const rowObj = [].reduce.call(ths, (acc, th, i) => {
            const { textContent } = row.cells[i];
            acc[th.textContent] = textContent;
            return acc;
        }, {});

        /* update row style */
        if (!row.classList.contains(selectedClass)) {
            var pastRow = document.getElementById('lora_summary').getElementsByClassName(selectedClass)[0];
            if (pastRow) {
                /* hide past row */
                pastRow.classList.remove(selectedClass);
                const pastRowObj = [].reduce.call(ths, (acc, th, i) => {
                    const { textContent } = pastRow.cells[i];
                    acc[th.textContent] = textContent;
                    return acc;
                }, {});
                document.getElementById(pastRowObj.Sample).style.display = "none";
            }
            row.classList.add(selectedClass);
            document.getElementById(rowObj.Sample).style.display = "block";
        }
    });
});

/** open the first contig information*/
document.getElementById("defaultOpenAssembly").click();

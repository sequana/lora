document.body.style.setProperty("--calc-svg-height", "auto");
const svg = document.getElementsByTagName('svg')[0];
document.body.style.setProperty("--calc-svg-height", svg.getBBox().height + "pt");

document.body.style.setProperty("--calc-dep-height", "auto");
const tableHeight = document.getElementById('dependencies').firstElementChild.clientHeight + "px";
document.body.style.setProperty("--calc-dep-height", tableHeight);

/* Hide div */
const hideDiv = (currentElm, idToHide) => {
  const elementToHide = document.getElementById(idToHide);

  if (!elementToHide.classList.contains(`open-section-${idToHide}`)) {
    elementToHide.classList.add(`open-section-${idToHide}`);
  } else {
    elementToHide.classList.remove(`open-section-${idToHide}`);
  }
  var spanElem = currentElm.querySelector('span');
  spanElem.textContent = spanElem.textContent === "expand_more" ? "expand_less" : "expand_more";
};

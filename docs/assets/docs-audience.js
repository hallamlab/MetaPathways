/* Reveal optional settings when following an existing section link. */
(() => {
  function revealLinkedSection() {
    if (!window.location.hash) return;
    let id;
    try { id = decodeURIComponent(window.location.hash.slice(1)); }
    catch (_) { return; }
    const target = document.getElementById(id);
    if (!target) return;
    let panel = target.closest('details');
    let opened = false;
    while (panel) {
      panel.open = true;
      opened = true;
      panel = panel.parentElement.closest('details');
    }
    if (opened) target.scrollIntoView();
  }
  document.addEventListener('DOMContentLoaded', revealLinkedSection);
  window.addEventListener('hashchange', revealLinkedSection);
  let printPanels = [];
  window.addEventListener('beforeprint', () => {
    printPanels = [...document.querySelectorAll('details.mp-advanced:not([open])')];
    printPanels.forEach(panel => { panel.open = true; });
  });
  window.addEventListener('afterprint', () => {
    printPanels.forEach(panel => { panel.open = false; });
    printPanels = [];
  });
})();

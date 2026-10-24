// Shared by every page: FAQ dialog, remembered chromosome, small helpers.
const SITE = (() => {
  const faq = document.getElementById("faq");
  document.getElementById("faqBtn").addEventListener("click", () => faq.showModal());
  document.getElementById("faqClose").addEventListener("click", () => faq.close());
  faq.addEventListener("click", e => { if (e.target === faq) faq.close(); });
  const store = {
    get(k, d) { try { const v = localStorage.getItem("idiograms:" + k); return v === null ? d : JSON.parse(v); } catch (e) { return d; } },
    set(k, v) { try { localStorage.setItem("idiograms:" + k, JSON.stringify(v)); } catch (e) {} },
  };
  const ORDER = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"];
  // chromosome from #chrN, else the one last opened on any page
  function initialChr(fallback) {
    const h = (location.hash || "").replace(/^#chr/i, "").toUpperCase();
    if (ORDER.includes(h)) return h;
    const s = store.get("chr", null);
    return ORDER.includes(s) ? s : fallback;
  }
  function rememberChr(c) { store.set("chr", c); try { history.replaceState(null, "", "#chr" + c); } catch (e) {} }
  // keep the remembered chromosome on links between tools
  document.querySelectorAll(".tabs a[data-chr-link]").forEach(a => a.addEventListener("click", () => {
    const c = store.get("chr", null); if (ORDER.includes(c) && !a.href.includes("#")) a.href += "#chr" + c;
  }));
  return { store, ORDER, initialChr, rememberChr };
})();

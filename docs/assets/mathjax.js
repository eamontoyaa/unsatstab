window.MathJax = {
  tex: {
    inlineMath: [["$", "$"], ["\\(", "\\)"]],
    displayMath: [["$$", "$$"], ["\\[", "\\]"]],
    processEscapes: true
  },
  options: {
    // Ignore explicitly marked elements, not every nested HTML element.
    ignoreHtmlClass: "mathjax_ignore|tex2jax_ignore",
    // Process both MkDocs math wrappers and notebook Markdown cells.
    processHtmlClass: "arithmatex|math|jp-RenderedMarkdown"
  }
};

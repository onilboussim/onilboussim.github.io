/**
 * Theme toggle — dark / light mode
 * Persists preference via localStorage
 * Falls back to system preference on first visit
 */
(function () {
  const STORAGE_KEY = 'ob-theme';

  // Determine initial theme: saved preference → system preference → light
  function getInitialTheme() {
    const saved = localStorage.getItem(STORAGE_KEY);
    if (saved === 'dark' || saved === 'light') return saved;
    if (window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches) {
      return 'dark';
    }
    return 'light';
  }

  // Apply theme to <html> element
  function applyTheme(theme) {
    document.documentElement.setAttribute('data-theme', theme);
    localStorage.setItem(STORAGE_KEY, theme);
  }

  // Apply immediately (before paint) to avoid flash
  applyTheme(getInitialTheme());

  // Wire up the toggle button after DOM is ready
  document.addEventListener('DOMContentLoaded', function () {
    var btn = document.getElementById('theme-toggle');
    if (!btn) return;

    btn.addEventListener('click', function () {
      var current = document.documentElement.getAttribute('data-theme');
      applyTheme(current === 'dark' ? 'light' : 'dark');
    });
  });
})();

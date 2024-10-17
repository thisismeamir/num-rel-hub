/*
THIS IS A GENERATED/BUNDLED FILE BY ESBUILD
if you want to view the source, please visit the github repository of this plugin
*/

var __defProp = Object.defineProperty;
var __getOwnPropDesc = Object.getOwnPropertyDescriptor;
var __getOwnPropNames = Object.getOwnPropertyNames;
var __hasOwnProp = Object.prototype.hasOwnProperty;
var __export = (target, all) => {
  for (var name in all)
    __defProp(target, name, { get: all[name], enumerable: true });
};
var __copyProps = (to, from, except, desc) => {
  if (from && typeof from === "object" || typeof from === "function") {
    for (let key of __getOwnPropNames(from))
      if (!__hasOwnProp.call(to, key) && key !== except)
        __defProp(to, key, { get: () => from[key], enumerable: !(desc = __getOwnPropDesc(from, key)) || desc.enumerable });
  }
  return to;
};
var __toCommonJS = (mod) => __copyProps(__defProp({}, "__esModule", { value: true }), mod);

// main.ts
var main_exports = {};
__export(main_exports, {
  default: () => ObsidianNiftyLinksPlugin,
  t: () => t
});
module.exports = __toCommonJS(main_exports);
var import_obsidian = require("obsidian");

// locale/en.ts
var en_default = {
  "Nifty Links": "Nifty Links",
  "Fixed width": "Fixed width",
  "Set the width of Nifty Links cards to a fixed 700px": "Set the width of Nifty Links cards to a fixed 700px",
  "Convert to Nifty Link": "Convert to Nifty Link"
};

// locale/zh-cn.ts
var zh_cn_default = {
  "Nifty Links": "Nifty Links",
  "Fixed width": "\u56FA\u5B9A\u5BBD\u5EA6",
  "Set the width of Nifty Links cards to a fixed 700px": "\u4FDD\u6301\u94FE\u63A5\u5361\u7247\u7684\u680F\u5BBD(700px)",
  "Convert to Nifty Link": "\u8F6C\u6362\u4E3ANifty Link"
};

// main.ts
var import_obsidian2 = require("obsidian");
var localeMap = {
  en: en_default,
  "zh-cn": zh_cn_default
};
var locale = localeMap[import_obsidian.moment.locale()];
function t(str) {
  return locale && locale[str] || en_default[str];
}
var DEFAULT_SETTINGS = {
  fixedWidth: false
};
var NiftyLinksSettingTab = class extends import_obsidian2.PluginSettingTab {
  constructor(app, plugin) {
    super(app, plugin);
    this.plugin = plugin;
  }
  async display() {
    const { containerEl } = this;
    containerEl.empty();
    new import_obsidian2.Setting(containerEl).setName(t("Fixed width")).setDesc(t("Set the width of Nifty Links cards to a fixed 700px")).addToggle((toggle) => toggle.setValue(this.plugin.settings.fixedWidth).onChange(async (value) => {
      this.plugin.settings.fixedWidth = value;
      await this.plugin.saveSettings();
      this.plugin.updateStyles();
    }));
  }
};
var ObsidianNiftyLinksPlugin = class extends import_obsidian.Plugin {
  async onload() {
    await this.loadSettings();
    this.addSettingTab(new NiftyLinksSettingTab(this.app, this));
    this.registerEvent(
      this.app.workspace.on("editor-menu", (menu, editor) => {
        const selection = editor.getSelection();
        if (selection && this.isUrl(selection.trim())) {
          menu.addItem((item) => {
            item.setTitle(t("Convert to Nifty Link")).setIcon("link").onClick(async () => {
              await this.urlToMarkdown(editor);
            });
          });
        }
      })
    );
    this.registerMarkdownCodeBlockProcessor("NiftyLinks", (source, el, ctx) => {
      const data = source.split("\n").reduce((acc, line) => {
        const [key, ...value] = line.split(": ");
        acc[key.trim()] = value.join(": ").trim();
        return acc;
      }, {});
      const url = data.url;
      let title = data.title || "";
      let description = data.description || "";
      const imageLink = data.image;
      const iconLink = data.favicon;
      title = title.replace(/\s{3,}/g, " ").trim();
      description = description.replace(/\s{3,}/g, " ").trim();
      const cardTextStyle = imageLink ? "" : ' style="width: 100%;"';
      const iconHTML = iconLink ? `<img class="nifty-link-icon" src="${iconLink}">` : "";
      const imageContainerHTML = imageLink ? `
		  <div class="nifty-link-image-container">
			<div class="nifty-link-image" style="background-image: url('${imageLink}')"></div>
		  </div>` : "";
      const html = `
		  <div class="nifty-link-card-container">
			<a class="nifty-link-card" href="${url}" target="_blank">
			  <div class="nifty-link-card-text"${cardTextStyle}>
				<div class="nifty-link-card-title">${title}</div>
				<div class="nifty-link-card-description">${description}</div>
				<div class="nifty-link-href">
				  ${iconHTML}${url}
				</div>
			  </div>
			  ${imageContainerHTML}
			</a>
		  </div>
		`;
      el.innerHTML = html;
    });
    this.updateStyles();
  }
  onunload() {
    console.log("unloading plugin");
  }
  isUrl(text) {
    const urlRegex = new RegExp("^(http:\\/\\/www\\.|https:\\/\\/www\\.|http:\\/\\/|https:\\/\\/)?[a-z0-9]+([\\-.]{1}[a-z0-9]+)*\\.[a-z]{2,5}(:[0-9]{1,5})?(\\/.*)?$");
    return urlRegex.test(text);
  }
  async urlToMarkdown(editor) {
    var _a, _b;
    let selectedText = editor.somethingSelected() ? editor.getSelection().trim() : false;
    if (selectedText && this.isUrl(selectedText)) {
      let url = selectedText;
      let api = "";
      const specialDomains = ["medium.com"];
      let isSpecialDomain = specialDomains.some((domain) => url.includes(domain));
      if (isSpecialDomain) {
        api = `https://api.microlink.io/?url=${url}`;
      } else {
        api = `http://iframely.server.crestify.com/iframely?url=${url}`;
      }
      try {
        let response = await (0, import_obsidian.requestUrl)({ url: api });
        let data = isSpecialDomain ? response.json.data : response.json;
        if (!isSpecialDomain && data.code === 403) {
          api = `https://api.microlink.io/?url=${url}`;
          response = await (0, import_obsidian.requestUrl)({ url: api });
          data = response.json.data;
          isSpecialDomain = true;
        }
        const imageLink = isSpecialDomain ? data.image ? data.image.url : "" : ((_a = data.links.find((value) => value.type.startsWith("image") && value.rel.includes("twitter"))) == null ? void 0 : _a.href) || "";
        const iconLink = isSpecialDomain ? data.logo ? data.logo.url : "" : ((_b = data.links.find((value) => value.type.startsWith("image") && value.rel.includes("icon"))) == null ? void 0 : _b.href) || "";
        let markdownLink = `
\`\`\`NiftyLinks
url: ${isSpecialDomain ? data.url || url : url}
title: ${isSpecialDomain ? data.title : data.meta.title || ""}
description: ${isSpecialDomain ? data.description : data.meta.description || ""}
favicon: ${iconLink}
${imageLink ? `image: ${imageLink}` : ""}
\`\`\`
`;
        editor.replaceSelection(markdownLink);
        return true;
      } catch (error) {
        return false;
      }
    } else {
      return false;
    }
  }
  async loadSettings() {
    this.settings = Object.assign({}, DEFAULT_SETTINGS, await this.loadData());
  }
  async saveSettings() {
    await this.saveData(this.settings);
  }
  updateStyles() {
    document.body.classList.toggle("nifty-links-fixed-width", this.settings.fixedWidth);
  }
};
//# sourceMappingURL=data:application/json;base64,ewogICJ2ZXJzaW9uIjogMywKICAic291cmNlcyI6IFsibWFpbi50cyIsICJsb2NhbGUvZW4udHMiLCAibG9jYWxlL3poLWNuLnRzIl0sCiAgInNvdXJjZXNDb250ZW50IjogWyJpbXBvcnQge1xyXG4gICAgQXBwLFxyXG4gICAgRWRpdG9yLFxyXG4gICAgTWFya2Rvd25WaWV3LFxyXG4gICAgUGx1Z2luLFxyXG4gICAgbW9tZW50LFxyXG4gICAgcmVxdWVzdFVybCxcclxufSBmcm9tIFwib2JzaWRpYW5cIjtcclxuXHJcbmltcG9ydCBlbiBmcm9tIFwiLi9sb2NhbGUvZW5cIjtcclxuaW1wb3J0IHpoQ04gZnJvbSBcIi4vbG9jYWxlL3poLWNuXCI7XHJcblxyXG5jb25zdCBsb2NhbGVNYXA6IHsgW2s6IHN0cmluZ106IFBhcnRpYWw8dHlwZW9mIGVuPiB9ID0ge1xyXG4gIGVuLFxyXG4gIFwiemgtY25cIjogemhDTixcclxufTtcclxuXHJcbmNvbnN0IGxvY2FsZSA9IGxvY2FsZU1hcFttb21lbnQubG9jYWxlKCldO1xyXG5cclxuZXhwb3J0IGZ1bmN0aW9uIHQoc3RyOiBrZXlvZiB0eXBlb2YgZW4pOiBzdHJpbmcge1xyXG4gIHJldHVybiAobG9jYWxlICYmIGxvY2FsZVtzdHJdKSB8fCBlbltzdHJdO1xyXG59XHJcblxyXG5pbnRlcmZhY2UgT2JzaWRpYW5OaWZ0eUxpbmtzUGx1Z2luU2V0dGluZ3Mge1xyXG4gICAgZml4ZWRXaWR0aDogYm9vbGVhbjtcclxufVxyXG5cclxuY29uc3QgREVGQVVMVF9TRVRUSU5HUzogT2JzaWRpYW5OaWZ0eUxpbmtzUGx1Z2luU2V0dGluZ3MgPSB7XHJcbiAgICBmaXhlZFdpZHRoOiBmYWxzZVxyXG59O1xyXG5cclxuaW1wb3J0IHsgUGx1Z2luU2V0dGluZ1RhYiwgU2V0dGluZyB9IGZyb20gXCJvYnNpZGlhblwiO1xyXG5cclxuY2xhc3MgTmlmdHlMaW5rc1NldHRpbmdUYWIgZXh0ZW5kcyBQbHVnaW5TZXR0aW5nVGFiIHtcclxuICAgIHBsdWdpbjogT2JzaWRpYW5OaWZ0eUxpbmtzUGx1Z2luO1xyXG5cclxuICAgIGNvbnN0cnVjdG9yKGFwcDogQXBwLCBwbHVnaW46IE9ic2lkaWFuTmlmdHlMaW5rc1BsdWdpbikge1xyXG4gICAgICAgIHN1cGVyKGFwcCwgcGx1Z2luKTtcclxuICAgICAgICB0aGlzLnBsdWdpbiA9IHBsdWdpbjtcclxuICAgIH1cclxuXHJcbiAgICBhc3luYyBkaXNwbGF5KCk6IFByb21pc2U8dm9pZD4ge1xyXG4gICAgICAgIGNvbnN0IHsgY29udGFpbmVyRWwgfSA9IHRoaXM7XHJcblxyXG4gICAgICAgIGNvbnRhaW5lckVsLmVtcHR5KCk7XHJcblxyXG4gICAgICAgIG5ldyBTZXR0aW5nKGNvbnRhaW5lckVsKVxyXG4gICAgICAgICAgICAuc2V0TmFtZSh0KCdGaXhlZCB3aWR0aCcpKVxyXG4gICAgICAgICAgICAuc2V0RGVzYyh0KCdTZXQgdGhlIHdpZHRoIG9mIE5pZnR5IExpbmtzIGNhcmRzIHRvIGEgZml4ZWQgNzAwcHgnKSlcclxuICAgICAgICAgICAgLmFkZFRvZ2dsZSh0b2dnbGUgPT4gdG9nZ2xlXHJcbiAgICAgICAgICAgICAgICAuc2V0VmFsdWUodGhpcy5wbHVnaW4uc2V0dGluZ3MuZml4ZWRXaWR0aClcclxuICAgICAgICAgICAgICAgIC5vbkNoYW5nZShhc3luYyAodmFsdWUpID0+IHtcclxuICAgICAgICAgICAgICAgICAgICB0aGlzLnBsdWdpbi5zZXR0aW5ncy5maXhlZFdpZHRoID0gdmFsdWU7XHJcbiAgICAgICAgICAgICAgICAgICAgYXdhaXQgdGhpcy5wbHVnaW4uc2F2ZVNldHRpbmdzKCk7XHJcbiAgICAgICAgICAgICAgICAgICAgdGhpcy5wbHVnaW4udXBkYXRlU3R5bGVzKCk7XHJcbiAgICAgICAgICAgICAgICB9KSk7XHJcbiAgICB9XHJcbn1cclxuXHJcbmV4cG9ydCBkZWZhdWx0IGNsYXNzIE9ic2lkaWFuTmlmdHlMaW5rc1BsdWdpbiBleHRlbmRzIFBsdWdpbiB7XHJcbiAgICBzZXR0aW5nczogT2JzaWRpYW5OaWZ0eUxpbmtzUGx1Z2luU2V0dGluZ3M7XHJcblxyXG5cdGFzeW5jIG9ubG9hZCgpIHtcclxuXHJcblx0XHRhd2FpdCB0aGlzLmxvYWRTZXR0aW5ncygpO1xyXG5cclxuXHRcdHRoaXMuYWRkU2V0dGluZ1RhYihuZXcgTmlmdHlMaW5rc1NldHRpbmdUYWIodGhpcy5hcHAsIHRoaXMpKTtcclxuXHJcblx0XHR0aGlzLnJlZ2lzdGVyRXZlbnQoXHJcblx0XHRcdHRoaXMuYXBwLndvcmtzcGFjZS5vbihcImVkaXRvci1tZW51XCIsIChtZW51LCBlZGl0b3I6IEVkaXRvcikgPT4ge1xyXG5cdFx0XHRcdGNvbnN0IHNlbGVjdGlvbiA9IGVkaXRvci5nZXRTZWxlY3Rpb24oKTtcclxuXHRcdFx0XHRpZiAoc2VsZWN0aW9uICYmIHRoaXMuaXNVcmwoc2VsZWN0aW9uLnRyaW0oKSkpIHtcclxuXHRcdFx0XHRcdG1lbnUuYWRkSXRlbSgoaXRlbSkgPT4ge1xyXG5cdFx0XHRcdFx0XHRpdGVtXHJcblx0XHRcdFx0XHRcdFx0LnNldFRpdGxlKHQoXCJDb252ZXJ0IHRvIE5pZnR5IExpbmtcIikpXHJcblx0XHRcdFx0XHRcdFx0LnNldEljb24oXCJsaW5rXCIpXHJcblx0XHRcdFx0XHRcdFx0Lm9uQ2xpY2soYXN5bmMgKCkgPT4ge1xyXG5cdFx0XHRcdFx0XHRcdFx0YXdhaXQgdGhpcy51cmxUb01hcmtkb3duKGVkaXRvcik7XHJcblx0XHRcdFx0XHRcdFx0fSk7XHJcblx0XHRcdFx0XHR9KTtcclxuXHRcdFx0XHR9XHJcblx0XHRcdH0pXHJcblx0XHQpO1xyXG5cclxuXHRcdHRoaXMucmVnaXN0ZXJNYXJrZG93bkNvZGVCbG9ja1Byb2Nlc3NvcihcIk5pZnR5TGlua3NcIiwgKHNvdXJjZSwgZWwsIGN0eCkgPT4ge1xyXG5cdFx0XHRjb25zdCBkYXRhID0gc291cmNlLnNwbGl0KCdcXG4nKS5yZWR1Y2UoKGFjYywgbGluZSkgPT4ge1xyXG5cdFx0XHRcdGNvbnN0IFtrZXksIC4uLnZhbHVlXSA9IGxpbmUuc3BsaXQoJzogJyk7XHJcblx0XHRcdFx0YWNjW2tleS50cmltKCldID0gdmFsdWUuam9pbignOiAnKS50cmltKCk7XHJcblx0XHRcdFx0cmV0dXJuIGFjYztcclxuXHRcdFx0fSwge30pO1xyXG5cclxuXHRcdFx0Y29uc3QgdXJsID0gZGF0YS51cmw7XHJcblx0XHRcdGxldCB0aXRsZSA9IGRhdGEudGl0bGUgfHwgXCJcIjtcclxuXHRcdFx0bGV0IGRlc2NyaXB0aW9uID0gZGF0YS5kZXNjcmlwdGlvbiB8fCBcIlwiO1xyXG5cdFx0XHRjb25zdCBpbWFnZUxpbmsgPSBkYXRhLmltYWdlO1xyXG5cdFx0XHRjb25zdCBpY29uTGluayA9IGRhdGEuZmF2aWNvbjtcclxuXHJcblx0XHRcdHRpdGxlID0gdGl0bGUucmVwbGFjZSgvXFxzezMsfS9nLCAnICcpLnRyaW0oKTtcclxuXHRcdFx0ZGVzY3JpcHRpb24gPSBkZXNjcmlwdGlvbi5yZXBsYWNlKC9cXHN7Myx9L2csICcgJykudHJpbSgpO1xyXG5cclxuXHRcdFx0Y29uc3QgY2FyZFRleHRTdHlsZSA9IGltYWdlTGluayA/IFwiXCIgOiAnIHN0eWxlPVwid2lkdGg6IDEwMCU7XCInO1xyXG5cclxuXHRcdFx0Y29uc3QgaWNvbkhUTUwgPSBpY29uTGluayA/IGA8aW1nIGNsYXNzPVwibmlmdHktbGluay1pY29uXCIgc3JjPVwiJHtpY29uTGlua31cIj5gIDogJyc7XHJcblxyXG5cdFx0XHRjb25zdCBpbWFnZUNvbnRhaW5lckhUTUwgPSBpbWFnZUxpbmsgPyBgXHJcblx0XHQgIDxkaXYgY2xhc3M9XCJuaWZ0eS1saW5rLWltYWdlLWNvbnRhaW5lclwiPlxyXG5cdFx0XHQ8ZGl2IGNsYXNzPVwibmlmdHktbGluay1pbWFnZVwiIHN0eWxlPVwiYmFja2dyb3VuZC1pbWFnZTogdXJsKCcke2ltYWdlTGlua30nKVwiPjwvZGl2PlxyXG5cdFx0ICA8L2Rpdj5gIDogJyc7XHJcblxyXG5cdFx0XHRjb25zdCBodG1sID0gYFxyXG5cdFx0ICA8ZGl2IGNsYXNzPVwibmlmdHktbGluay1jYXJkLWNvbnRhaW5lclwiPlxyXG5cdFx0XHQ8YSBjbGFzcz1cIm5pZnR5LWxpbmstY2FyZFwiIGhyZWY9XCIke3VybH1cIiB0YXJnZXQ9XCJfYmxhbmtcIj5cclxuXHRcdFx0ICA8ZGl2IGNsYXNzPVwibmlmdHktbGluay1jYXJkLXRleHRcIiR7Y2FyZFRleHRTdHlsZX0+XHJcblx0XHRcdFx0PGRpdiBjbGFzcz1cIm5pZnR5LWxpbmstY2FyZC10aXRsZVwiPiR7dGl0bGV9PC9kaXY+XHJcblx0XHRcdFx0PGRpdiBjbGFzcz1cIm5pZnR5LWxpbmstY2FyZC1kZXNjcmlwdGlvblwiPiR7ZGVzY3JpcHRpb259PC9kaXY+XHJcblx0XHRcdFx0PGRpdiBjbGFzcz1cIm5pZnR5LWxpbmstaHJlZlwiPlxyXG5cdFx0XHRcdCAgJHtpY29uSFRNTH0ke3VybH1cclxuXHRcdFx0XHQ8L2Rpdj5cclxuXHRcdFx0ICA8L2Rpdj5cclxuXHRcdFx0ICAke2ltYWdlQ29udGFpbmVySFRNTH1cclxuXHRcdFx0PC9hPlxyXG5cdFx0ICA8L2Rpdj5cclxuXHRcdGA7XHJcblxyXG5cdFx0XHRlbC5pbm5lckhUTUwgPSBodG1sO1xyXG5cdFx0fSk7XHJcblxyXG5cdFx0dGhpcy51cGRhdGVTdHlsZXMoKTtcclxuXHR9XHJcblxyXG5cdG9udW5sb2FkKCkge1xyXG5cdFx0Y29uc29sZS5sb2coXCJ1bmxvYWRpbmcgcGx1Z2luXCIpO1xyXG5cdH1cclxuXHJcblx0aXNVcmwodGV4dCkge1xyXG5cdFx0Y29uc3QgdXJsUmVnZXggPSBuZXcgUmVnRXhwKFwiXihodHRwOlxcXFwvXFxcXC93d3dcXFxcLnxodHRwczpcXFxcL1xcXFwvd3d3XFxcXC58aHR0cDpcXFxcL1xcXFwvfGh0dHBzOlxcXFwvXFxcXC8pP1thLXowLTldKyhbXFxcXC0uXXsxfVthLXowLTldKykqXFxcXC5bYS16XXsyLDV9KDpbMC05XXsxLDV9KT8oXFxcXC8uKik/JFwiKTtcclxuXHRcdHJldHVybiB1cmxSZWdleC50ZXN0KHRleHQpO1xyXG5cdH1cclxuXHJcblx0YXN5bmMgdXJsVG9NYXJrZG93bihlZGl0b3IpIHtcclxuXHRcdGxldCBzZWxlY3RlZFRleHQgPSBlZGl0b3Iuc29tZXRoaW5nU2VsZWN0ZWQoKVxyXG5cdFx0XHQ/IGVkaXRvci5nZXRTZWxlY3Rpb24oKS50cmltKClcclxuXHRcdFx0OiBmYWxzZTtcclxuXHRcdGlmIChzZWxlY3RlZFRleHQgJiYgdGhpcy5pc1VybChzZWxlY3RlZFRleHQpKSB7XHJcblx0XHRcdGxldCB1cmwgPSBzZWxlY3RlZFRleHQ7XHJcblx0XHRcdGxldCBhcGkgPSBcIlwiO1xyXG5cdFx0XHRjb25zdCBzcGVjaWFsRG9tYWlucyA9IFtcIm1lZGl1bS5jb21cIl07XHJcblx0XHRcdGxldCBpc1NwZWNpYWxEb21haW4gPSBzcGVjaWFsRG9tYWlucy5zb21lKGRvbWFpbiA9PiB1cmwuaW5jbHVkZXMoZG9tYWluKSk7XHJcblx0XHRcdGlmIChpc1NwZWNpYWxEb21haW4pIHtcclxuXHRcdFx0XHRhcGkgPSBgaHR0cHM6Ly9hcGkubWljcm9saW5rLmlvLz91cmw9JHt1cmx9YDtcclxuXHRcdFx0fSBlbHNlIHtcclxuXHRcdFx0XHRhcGkgPSBgaHR0cDovL2lmcmFtZWx5LnNlcnZlci5jcmVzdGlmeS5jb20vaWZyYW1lbHk/dXJsPSR7dXJsfWA7XHJcblx0XHRcdH1cclxuXHJcblx0XHRcdHRyeSB7XHJcblx0XHRcdFx0bGV0IHJlc3BvbnNlID0gYXdhaXQgcmVxdWVzdFVybCh7IHVybDogYXBpIH0pO1xyXG5cdFx0XHRcdGxldCBkYXRhID0gaXNTcGVjaWFsRG9tYWluID8gcmVzcG9uc2UuanNvbi5kYXRhIDogcmVzcG9uc2UuanNvbjtcclxuXHRcdFx0XHRpZiAoIWlzU3BlY2lhbERvbWFpbiAmJiBkYXRhLmNvZGUgPT09IDQwMykge1xyXG5cdFx0XHRcdFx0YXBpID0gYGh0dHBzOi8vYXBpLm1pY3JvbGluay5pby8/dXJsPSR7dXJsfWA7XHJcblx0XHRcdFx0XHRyZXNwb25zZSA9IGF3YWl0IHJlcXVlc3RVcmwoeyB1cmw6IGFwaSB9KTtcclxuXHRcdFx0XHRcdGRhdGEgPSByZXNwb25zZS5qc29uLmRhdGE7XHJcblx0XHRcdFx0XHRpc1NwZWNpYWxEb21haW4gPSB0cnVlO1xyXG5cdFx0XHRcdH1cclxuXHRcdFx0XHRjb25zdCBpbWFnZUxpbmsgPSBpc1NwZWNpYWxEb21haW4gPyAoZGF0YS5pbWFnZSA/IGRhdGEuaW1hZ2UudXJsIDogJycpIDogZGF0YS5saW5rcy5maW5kKCh2YWx1ZSkgPT4gdmFsdWUudHlwZS5zdGFydHNXaXRoKFwiaW1hZ2VcIikgJiYgdmFsdWUucmVsLmluY2x1ZGVzKCd0d2l0dGVyJykpPy5ocmVmIHx8ICcnO1xyXG5cdFx0XHRcdGNvbnN0IGljb25MaW5rID0gaXNTcGVjaWFsRG9tYWluID8gKGRhdGEubG9nbyA/IGRhdGEubG9nby51cmwgOiAnJykgOiBkYXRhLmxpbmtzLmZpbmQoKHZhbHVlKSA9PiB2YWx1ZS50eXBlLnN0YXJ0c1dpdGgoXCJpbWFnZVwiKSAmJiB2YWx1ZS5yZWwuaW5jbHVkZXMoJ2ljb24nKSk/LmhyZWYgfHwgJyc7XHJcblxyXG5cdFx0XHRcdGxldCBtYXJrZG93bkxpbmsgPSBgXFxuXFxgXFxgXFxgTmlmdHlMaW5rc1xyXG51cmw6ICR7aXNTcGVjaWFsRG9tYWluID8gKGRhdGEudXJsIHx8IHVybCkgOiB1cmx9XHJcbnRpdGxlOiAke2lzU3BlY2lhbERvbWFpbiA/IGRhdGEudGl0bGUgOiBkYXRhLm1ldGEudGl0bGUgfHwgXCJcIn1cclxuZGVzY3JpcHRpb246ICR7aXNTcGVjaWFsRG9tYWluID8gZGF0YS5kZXNjcmlwdGlvbiA6IGRhdGEubWV0YS5kZXNjcmlwdGlvbiB8fCBcIlwifVxyXG5mYXZpY29uOiAke2ljb25MaW5rfVxyXG4ke2ltYWdlTGluayA/IGBpbWFnZTogJHtpbWFnZUxpbmt9YCA6IFwiXCJ9XHJcblxcYFxcYFxcYFxcbmA7XHJcblxyXG4gICAgICAgICAgICBlZGl0b3IucmVwbGFjZVNlbGVjdGlvbihtYXJrZG93bkxpbmspO1xyXG4gICAgICAgICAgICByZXR1cm4gdHJ1ZTtcclxuICAgICAgICB9IGNhdGNoIChlcnJvcikge1xyXG4gICAgICAgICAgICByZXR1cm4gZmFsc2U7XHJcbiAgICAgICAgfVxyXG4gICAgfSBlbHNlIHtcclxuICAgICAgICByZXR1cm4gZmFsc2U7XHJcbiAgICB9XHJcblx0fVxyXG5cclxuXHRhc3luYyBsb2FkU2V0dGluZ3MoKSB7XHJcblx0XHR0aGlzLnNldHRpbmdzID0gT2JqZWN0LmFzc2lnbih7fSwgREVGQVVMVF9TRVRUSU5HUywgYXdhaXQgdGhpcy5sb2FkRGF0YSgpKTtcclxuXHR9XHJcblxyXG5cdGFzeW5jIHNhdmVTZXR0aW5ncygpIHtcclxuXHRcdGF3YWl0IHRoaXMuc2F2ZURhdGEodGhpcy5zZXR0aW5ncyk7XHJcblx0fVxyXG5cclxuXHR1cGRhdGVTdHlsZXMoKSB7XHJcblx0XHRkb2N1bWVudC5ib2R5LmNsYXNzTGlzdC50b2dnbGUoJ25pZnR5LWxpbmtzLWZpeGVkLXdpZHRoJywgdGhpcy5zZXR0aW5ncy5maXhlZFdpZHRoKTtcclxuXHR9XHJcbn1cclxuIiwgImV4cG9ydCBkZWZhdWx0IHsgICAgXHJcbiAgICBcIk5pZnR5IExpbmtzXCI6IFwiTmlmdHkgTGlua3NcIixcclxuICAgIFwiRml4ZWQgd2lkdGhcIjogXCJGaXhlZCB3aWR0aFwiLFxyXG4gICAgXCJTZXQgdGhlIHdpZHRoIG9mIE5pZnR5IExpbmtzIGNhcmRzIHRvIGEgZml4ZWQgNzAwcHhcIjogXCJTZXQgdGhlIHdpZHRoIG9mIE5pZnR5IExpbmtzIGNhcmRzIHRvIGEgZml4ZWQgNzAwcHhcIixcclxuICAgIFwiQ29udmVydCB0byBOaWZ0eSBMaW5rXCI6IFwiQ29udmVydCB0byBOaWZ0eSBMaW5rXCJcclxufTsiLCAiZXhwb3J0IGRlZmF1bHQge1xyXG4gICAgXCJOaWZ0eSBMaW5rc1wiOiBcIk5pZnR5IExpbmtzXCIsXHJcbiAgICBcIkZpeGVkIHdpZHRoXCI6IFwiXHU1NkZBXHU1QjlBXHU1QkJEXHU1RUE2XCIsXHJcbiAgICBcIlNldCB0aGUgd2lkdGggb2YgTmlmdHkgTGlua3MgY2FyZHMgdG8gYSBmaXhlZCA3MDBweFwiOiBcIlx1NEZERFx1NjMwMVx1OTRGRVx1NjNBNVx1NTM2MVx1NzI0N1x1NzY4NFx1NjgwRlx1NUJCRCg3MDBweClcIixcclxuICAgIFwiQ29udmVydCB0byBOaWZ0eSBMaW5rXCI6IFwiXHU4RjZDXHU2MzYyXHU0RTNBTmlmdHkgTGlua1wiXHJcbn07Il0sCiAgIm1hcHBpbmdzIjogIjs7Ozs7Ozs7Ozs7Ozs7Ozs7Ozs7Ozs7O0FBQUE7QUFBQTtBQUFBO0FBQUE7QUFBQTtBQUFBO0FBQUEsc0JBT087OztBQ1BQLElBQU8sYUFBUTtBQUFBLEVBQ1gsZUFBZTtBQUFBLEVBQ2YsZUFBZTtBQUFBLEVBQ2YsdURBQXVEO0FBQUEsRUFDdkQseUJBQXlCO0FBQzdCOzs7QUNMQSxJQUFPLGdCQUFRO0FBQUEsRUFDWCxlQUFlO0FBQUEsRUFDZixlQUFlO0FBQUEsRUFDZix1REFBdUQ7QUFBQSxFQUN2RCx5QkFBeUI7QUFDN0I7OztBRjBCQSxJQUFBQSxtQkFBMEM7QUFuQjFDLElBQU0sWUFBaUQ7QUFBQSxFQUNyRDtBQUFBLEVBQ0EsU0FBUztBQUNYO0FBRUEsSUFBTSxTQUFTLFVBQVUsdUJBQU8sT0FBTyxDQUFDO0FBRWpDLFNBQVMsRUFBRSxLQUE4QjtBQUM5QyxTQUFRLFVBQVUsT0FBTyxHQUFHLEtBQU0sV0FBRyxHQUFHO0FBQzFDO0FBTUEsSUFBTSxtQkFBcUQ7QUFBQSxFQUN2RCxZQUFZO0FBQ2hCO0FBSUEsSUFBTSx1QkFBTixjQUFtQyxrQ0FBaUI7QUFBQSxFQUdoRCxZQUFZLEtBQVUsUUFBa0M7QUFDcEQsVUFBTSxLQUFLLE1BQU07QUFDakIsU0FBSyxTQUFTO0FBQUEsRUFDbEI7QUFBQSxFQUVBLE1BQU0sVUFBeUI7QUFDM0IsVUFBTSxFQUFFLFlBQVksSUFBSTtBQUV4QixnQkFBWSxNQUFNO0FBRWxCLFFBQUkseUJBQVEsV0FBVyxFQUNsQixRQUFRLEVBQUUsYUFBYSxDQUFDLEVBQ3hC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QUFNLG9EQUFvRDtBQUFBLE1BQzNEO0FBRUEsVUFBSTtBQUNILFlBQUksV0FBVyxVQUFNLDRCQUFXLEVBQUUsS0FBSyxJQUFJLENBQUM7QUFDNUMsWUFBSSxPQUFPLGtCQUFrQixTQUFTLEtBQUssT0FBTyxTQUFTO0FBQzNELFlBQUksQ0FBQyxtQkFBbUIsS0FBSyxTQUFTLEtBQUs7QUFDMUMsZ0JBQU0saUNBQWlDO0FBQ3ZDLHFCQUFXLFVBQU0sNEJBQVcsRUFBRSxLQUFLLElBQUksQ0FBQztBQUN4QyxpQkFBTyxTQUFTLEtBQUs7QUFDckIsNEJBQWtCO0FBQUEsUUFDbkI7QUFDQSxjQUFNLFlBQVksa0JBQW1CLEtBQUssUUFBUSxLQUFLLE1BQU0sTUFBTSxPQUFNLFVBQUssTUFBTSxLQUFLLENBQUMsVUFBVSxNQUFNLEtBQUssV0FBVyxPQUFPLEtBQUssTUFBTSxJQUFJLFNBQVMsU0FBUyxDQUFDLE1BQTFGLG1CQUE2RixTQUFRO0FBQzlLLGNBQU0sV0FBVyxrQkFBbUIsS0FBSyxPQUFPLEtBQUssS0FBSyxNQUFNLE9BQU0sVUFBSyxNQUFNLEtBQUssQ0FBQyxVQUFVLE1BQU0sS0FBSyxXQUFXLE9BQU8sS0FBSyxNQUFNLElBQUksU0FBUyxNQUFNLENBQUMsTUFBdkYsbUJBQTBGLFNBQVE7QUFFeEssWUFBSSxlQUFlO0FBQUE7QUFBQSxPQUNoQixrQkFBbUIsS0FBSyxPQUFPLE1BQU87QUFBQSxTQUNwQyxrQkFBa0IsS0FBSyxRQUFRLEtBQUssS0FBSyxTQUFTO0FBQUEsZUFDNUMsa0JBQWtCLEtBQUssY0FBYyxLQUFLLEtBQUssZUFBZTtBQUFBLFdBQ2xFO0FBQUEsRUFDVCxZQUFZLFVBQVUsY0FBYztBQUFBO0FBQUE7QUFHMUIsZUFBTyxpQkFBaUIsWUFBWTtBQUNwQyxlQUFPO0FBQUEsTUFDWCxTQUFTLE9BQVA7QUFDRSxlQUFPO0FBQUEsTUFDWDtBQUFBLElBQ0osT0FBTztBQUNILGFBQU87QUFBQSxJQUNYO0FBQUEsRUFDSDtBQUFBLEVBRUEsTUFBTSxlQUFlO0FBQ3BCLFNBQUssV0FBVyxPQUFPLE9BQU8sQ0FBQyxHQUFHLGtCQUFrQixNQUFNLEtBQUssU0FBUyxDQUFDO0FBQUEsRUFDMUU7QUFBQSxFQUVBLE1BQU0sZUFBZTtBQUNwQixVQUFNLEtBQUssU0FBUyxLQUFLLFFBQVE7QUFBQSxFQUNsQztBQUFBLEVBRUEsZUFBZTtBQUNkLGFBQVMsS0FBSyxVQUFVLE9BQU8sMkJBQTJCLEtBQUssU0FBUyxVQUFVO0FBQUEsRUFDbkY7QUFDRDsiLAogICJuYW1lcyI6IFsiaW1wb3J0X29ic2lkaWFuIl0KfQo=
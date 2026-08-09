# Minimal Scientific Web App Style Guide

This document is an implementation specification for Codex and other coding agents. Use it when creating or restyling small scientific web applications in the Edwards Lab ecosystem. The goal is a modern, understated interface that keeps the scientific task—not the visual design—at the center.

## Design intent

Build interfaces that feel calm, trustworthy, compact, and immediately usable. Prefer familiar application patterns over marketing-page decoration. The interface should resemble a well-made scientific utility: clear hierarchy, restrained color, visible state, and no unnecessary ornament.

Use these principles in priority order:

1. Make the primary task obvious within one screen.
2. Preserve user trust with explicit privacy and processing-state messages.
3. Keep controls compact while maintaining comfortable touch targets.
4. Use color to communicate action or status, never as decoration.
5. Maintain accessible focus, contrast, labels, and reduced-motion behavior.

## Visual foundation

Use the operating system UI font stack. Do not load a web font unless a product requirement explicitly calls for one.

```css
:root {
  --background: #f8fafc;
  --surface: #ffffff;
  --surface-subtle: #f1f5f9;
  --border: #e2e8f0;
  --border-strong: #cbd5e1;
  --text: #1e293b;
  --muted: #64748b;
  --primary: #2563eb;
  --primary-hover: #1d4ed8;
  --primary-soft: #eff6ff;
  --success: #15803d;
  --success-soft: #f0fdf4;
  --warning: #d97706;
  --danger: #dc2626;
  --radius: 6px;
  --shadow: 0 1px 3px rgb(15 23 42 / 8%);
  font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
}
```

The page canvas is pale slate, content surfaces are white, and primary actions are blue. Use green for successful or private/local states, amber for warnings, and red for failures or destructive actions. Do not introduce extra brand colors without a functional reason.

## Layout

- Limit the main content area to approximately `1040px` or `1100px`.
- Use `16px` horizontal page padding on small screens and `32px` where space permits.
- Use a white header with a subtle bottom border and one-pixel shadow.
- Put each major workflow stage in a white card with a `1px` border, `6px` radius, and subtle shadow.
- Use a spacing scale based on `4px`: 4, 8, 12, 16, 24, 32, and 48px.
- Keep introductory content short. A utility should not begin with a full-screen hero.
- Collapse multi-column controls to one column below approximately `720px`.

## Typography

- Body text: `14px` to `16px`, line-height `1.5` to `1.6`.
- Page title: responsive `32px` to `48px`, weight `700`, tight letter spacing.
- Card or section heading: approximately `17px`, weight `600` or `650`.
- Labels and buttons: `13px` to `14px`, weight `600`.
- Supporting copy: `12px` to `14px`, using the muted text color.
- File names, extensions, identifiers, and logs: system monospace.
- Do not mix serif and sans-serif display faces in this style.

## Components

### Header

Use a compact white header. Place the product name at the left and no more than two utility items at the right. A local-processing or privacy state may appear as a small green pill. Avoid oversized logos and complex navigation.

### Cards

Cards use a white surface, fine slate border, `6px` radius, and minimal shadow. Avoid floating cards with large radii or deep shadows. Nested options should usually use borders and subtle fills rather than additional shadows.

### Buttons

The primary button uses solid blue with white text. Secondary buttons use a slate fill or white surface with a border. Buttons use a `6px` radius and concise labels. Do not use gradients. Disable unavailable actions visibly and retain at least `44px` effective touch height on mobile.

### Upload areas

Use a dashed slate border on a pale background. On hover or drag-over, change the border to primary blue and the background to pale blue. State the supported inputs and offer click-to-browse as well as drag-and-drop.

### Selection controls

Display complex format choices as bordered option tiles. A selected tile uses a pale blue background and blue border. Keep titles, descriptions, and extensions aligned consistently. Native checkboxes may be visually hidden only when the replacement retains keyboard focus and a clear selected state.

### Status and privacy messages

Use lightly tinted panels with matching borders. Privacy/local-processing uses pale green or blue; warnings use pale amber; errors use pale red. Keep text direct and specific. Never imply that data stays local unless the implementation truly performs no upload.

### Download rows and data lists

Use compact rows on a pale slate surface. Align the primary action to the right. Truncate long file names visually without altering the actual downloadable name. Use monospace for filenames and identifiers.

## Interaction rules

- Provide visible `:focus-visible` outlines of at least `3px`.
- Do not rely on hover to expose essential information or actions.
- Keep animation to color and border transitions of roughly `150ms`.
- Respect `prefers-reduced-motion` and disable smooth scrolling and nonessential transitions.
- Announce async state changes through an `aria-live` status region.
- Use semantic HTML: `button` for actions, `a` for navigation/downloads, `fieldset` and `legend` for grouped options.
- Preserve working element IDs and event hooks during visual-only redesigns.

## Content style

Use plain, task-oriented language. Prefer “Choose a GenBank file” over promotional copy. Tell users what happens, what is accepted, and what remains local. Button labels should begin with a verb: “Choose file,” “Generate files,” or “Download.” Avoid exclamation marks, vague reassurance, and decorative technical jargon.

## Codex implementation checklist

Before completing a web app that follows this guide, Codex must verify:

- the primary workflow is visible and understandable without documentation;
- the layout works at desktop and mobile widths;
- every interactive control has a keyboard-visible focus state;
- color contrast remains readable across normal, muted, selected, disabled, success, and error states;
- loading, empty, success, and error states are present where relevant;
- privacy claims match actual network and storage behavior;
- reduced-motion preferences are honored;
- existing behavioral tests still pass;
- the application introduces no external font, icon, or runtime dependency solely for appearance.

When an established application has stronger repository-specific design conventions, preserve those conventions and use this guide only to fill gaps.

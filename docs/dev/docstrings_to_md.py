import importlib
import inspect
import os
import pkgutil
import sys



def format_long_signature(sig_str, max_len=60, indent_spaces=4):
    """
    Auto line-breaks long argument lists using HTML for Markdown rendering.
    """

    if len(sig_str) <= max_len:
        return sig_str

    # Extract what is inside the parentheses
    if not (sig_str.startswith("(") and sig_str.endswith(")")):
        return sig_str

    content = sig_str[1:-1]
    # Split arguments handling simple parameter lists safely
    args = [arg.strip() for arg in content.split(",") if arg.strip()]

    if not args:
        return "()"

    # Build a multi-line HTML string compatible with Markdown headings
    indent = "&nbsp;" * indent_spaces
    formatted_args = f",<br>{indent}".join(args)
    return f"(<br>{indent}{formatted_args})"



def generate_markdown(package_name, output_file="API.md"):
    try:
        package = importlib.import_module(package_name)
    except ImportError as e:
        print(f"Error: Could not import pkg '{package_name}'. Details: {e}")
        return

    markdown_lines = []
    markdown_lines.append(f"# {package_name}\n")

    pkg_doc = inspect.getdoc(package)
    if pkg_doc:
        markdown_lines.append(f"\n```text\n{pkg_doc}\n```\n")
    markdown_lines.append("\n---\n")

    submodules = []
    if hasattr(package, "__path__"):
        pkg_walk = pkgutil.walk_packages(
            package.__path__, package.__name__ + ".")

        for _, mod_name, _ in pkg_walk:
            submodules.append(mod_name)
    else:
        submodules.append(package_name)

    for mod_name in sorted(submodules):
        try:
            mod = importlib.import_module(mod_name)
        except Exception as e:
            print(f"Warning: Skipping {mod_name} due to import error: {e}")
            continue

        markdown_lines.append(f"## {mod_name}\n")
        mod_doc = inspect.getdoc(mod)
        if mod_doc:
            markdown_lines.append(f"\n```text\n{mod_doc}\n```\n")
        markdown_lines.append("\n")

        members = inspect.getmembers(mod)
        has_members = False

        for name, obj in sorted(members):
            if name.startswith("_"):
                continue

            obj_module = getattr(obj, "__module__", None)
            if obj_module != mod_name:
                continue

            is_func = inspect.isfunction(obj) or inspect.isbuiltin(obj)
            is_class = inspect.isclass(obj)
            is_constant = (
                not is_func and not is_class and not inspect.ismodule(obj)
            )

            if is_func or is_class or is_constant:
                has_members = True

                if is_class:
                    m_type = " `[Class]`"
                    try:
                        sig = str(inspect.signature(obj.__init__))
                        sig = sig.replace(
                            "(self, ", "(").replace("(self)", "()")
                        sig = format_long_signature(sig)
                    except Exception:
                        sig = ""
                elif is_func:
                    m_type = " `[Function]`"
                    try:
                        sig = str(inspect.signature(obj))
                        sig = format_long_signature(sig)
                    except Exception:
                        sig = "()"
                else:
                    m_type = " `[Constant]`"
                    sig = ""

                # Write Member Heading with HTML-formatted signature wrapping
                markdown_lines.append(f"### {name}{sig}{m_type}\n")

                doc = inspect.getdoc(obj)
                if doc:
                    markdown_lines.append(f"\n```text\n{doc}\n```\n\n")
                else:
                    # no doc
                    markdown_lines.append("\n\n")

                if is_class:
                    # Get all functions and methods defined directly
                    # inside this class
                    class_members = inspect.getmembers(
                        obj,
                        predicate=lambda m: (
                            inspect.isfunction(m)
                            or inspect.ismethod(m)
                        )
                    )

                    for method_name, method_obj in sorted(class_members):
                        if method_name.startswith("_"):
                            continue

                        try:
                            m_sig = str(inspect.signature(method_obj))
                            m_sig = m_sig.replace(
                                "(self, ", "(").replace("(self)", "()")
                            m_sig = format_long_signature(
                                m_sig, indent_spaces=8)
                        except Exception:
                            m_sig = "()"

                        # Format sub-methods neatly under the class header
                        markdown_lines.append(
                            f"#### `{name}.{method_name}{m_sig}` `[Method]`\n")

                        method_doc = inspect.getdoc(method_obj)
                        if method_doc:
                            markdown_lines.append(
                                f"```text\n{method_doc}\n```\n\n")
                        else:
                            markdown_lines.append("\n\n")

        if not has_members:
            # no public member
            markdown_lines.append("\n\n")

        markdown_lines.append("---\n\n")

    with open(output_file, "w", encoding="utf-8") as f:
        f.writelines(markdown_lines)

    print(f"Generated API Markdown: {output_file}")



if __name__ == "__main__":
    generate_markdown("rawk")
